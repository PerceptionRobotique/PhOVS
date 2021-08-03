/*!
  \class CFeatureDefocusedLuminance.cpp
  \brief Class that defines the defocused image brightness visual feature

  for more details see
  G. Caron, RAL/ICRA 2021 submission
*/

#include <visp/vpMatrix.h>
#include <visp/vpHomogeneousMatrix.h>
#include <visp/vpDisplay.h>
#include <visp/vpPixelMeterConversion.h>
#include <visp/vpMeterPixelConversion.h>
#include <visp/vpImageConvert.h>

#include <PhOVS/CImageFilter.h>
#include <PhOVS/CFeatureDefocusedLuminance.h>

/*!
  Initialize the memory space requested for vpFeatureLuminance visual feature.
*/
void
CFeatureDefocusedLuminance::init()
{
  recompute_xy = false;
  nbParameters = 1;
  dim_s = 0 ;
  flags = NULL;
  pixInfo = NULL;

  //imWidth = imHeight = 0;

  if (flags == NULL)
    flags = new bool[nbParameters];
  for (int i = 0; i < nbParameters; i++) flags[i] = false;

  //default value rho (1 meters)
  Z = 1.;

  firstTimeIn =0 ;

  nbDOF = 6;
  for(int i = 0 ; i < 6 ; i++) DOF[i] = true;
}

void
CFeatureDefocusedLuminance::init(int _imHeight, int _imWidth, int _di, int _dj, int _nbri, int _nbrj, int _pas, vpImage<unsigned char> *Imask, double _Z, int _derivativeMaskSize)
{
  init() ;

  imWidth = _imWidth;
  imHeight = _imHeight;
	
  di = _di;
  dj = _dj;
  if((_nbri == di) && (imHeight > di))
    nbri = imHeight - di;
  else
    nbri = _nbri;
  if((_nbrj == dj) && (imWidth > dj))
    nbrj = imWidth - dj;
  else
    nbrj = _nbrj;

  pas =_pas;

  // number of feature = nb column x nb lines in the images
  dim_s = (nbri/pas)*(nbrj/pas);

  if(Imask != NULL)
  {
    unsigned char *pt_Imask = Imask->bitmap;
    long int i, j, nbi = imHeight, nbj = imWidth;

    for(i = 0 ; i < di ; i+=pas)
      for(j = 0 ; j < nbj ; j+=pas, pt_Imask+=pas);
        
    for(i = di ; i < (di+nbri) ; i+=pas)
    {
      for(j = 0 ; j < dj ; j+=pas, pt_Imask+=pas);

      for(j = dj ; j < (dj+nbrj) ; j+=pas, pt_Imask+=pas)
        if((*pt_Imask)==0)
          dim_s--;

      for(j = dj+nbrj ; j < nbj ; j+=pas, pt_Imask+=pas);
    }
  }

  std::cout << "dim_s : " << dim_s << std::endl;
  s.resize(dim_s) ;
  
  if (pixInfo != NULL)
    delete [] pixInfo;

  nbDim = 2; // ?

  pixInfo = new CLuminanceDefocusCIP[dim_s];


  if(Imask != NULL)
  {
    unsigned char *pt_Imask = Imask->bitmap;
    long int i, j, nbi = imHeight, nbj = imWidth;
    
    for(i = 0 ; i < di ; i+=pas)
      for(j = 0 ; j < nbj ; j+=pas, pt_Imask+=pas);


        CLuminanceDefocusCIP *pt_pixInfo = (CLuminanceDefocusCIP *)pixInfo; 
        for(i = di ; i < (di+nbri) ; i+=pas)
        {
          for(j = 0 ; j < dj ; j+=pas, pt_Imask+=pas);

          for(j = dj ; j < (dj+nbrj) ; j+=pas, pt_Imask+=pas)
            if((*pt_Imask)!=0)
            {
              pt_pixInfo->i = i;
              pt_pixInfo->j = j;
              pt_pixInfo++;
            }

          for(j = (dj+nbrj) ; j < nbj ; j+=pas, pt_Imask+=pas);
        }
  }
  else
  {
    int i, j;


        CLuminanceDefocusCIP *pt_pixInfo = (CLuminanceDefocusCIP *)pixInfo; 
        for(i = di ; i < (di+nbri) ; i+=pas)
          for(j = dj ; j < (dj+nbrj) ; j+=pas)
          {
            pt_pixInfo->i = i;
            pt_pixInfo->j = j;
            pt_pixInfo++;
          }
  }
  
  Z = _Z ;

  //la taille est constante à 7 pour l'instant (7-1)/2 = 3
  derivativeMaskHalfSize = 3;//(_derivativeMaskSize-1)*0.5;
  nbNeigh = 2*derivativeMaskHalfSize;
}

/*! 
  Default constructor that build a visual feature.
*/
CFeatureDefocusedLuminance::CFeatureDefocusedLuminance() : vpBasicFeature(), di(10), dj(10), nbri(10), nbrj(10), pas(1), imWidth(0), imHeight(0), derivativeMaskHalfSize(3), nbNeigh(6), nbDim(0), cam(8e-3,5e-6,1,0.25,320,240), recompute_xy(false)
{
  init() ;
}

/*! 
  Default destructor.
*/
CFeatureDefocusedLuminance::~CFeatureDefocusedLuminance() 
{
  if (pixInfo != NULL)
  {
    delete [] pixInfo ; 
    pixInfo == NULL;
  }
  if (flags != NULL)
  {
    delete [] flags;
    flags = NULL;
  }
}


void
CFeatureDefocusedLuminance::set_DOF(bool un, bool deux, bool trois, bool quatre, bool cinq, bool six)
{
  nbDOF = 0;
  DOF[0] = un; if(un) nbDOF++;
  DOF[1] = deux;  if(deux) nbDOF++;
  DOF[2] = trois;  if(trois) nbDOF++;
  DOF[3] = quatre;  if(quatre) nbDOF++;
  DOF[4] = cinq;  if(cinq) nbDOF++;
  DOF[5] = six;  if(six) nbDOF++;
}

/*!
  Set the value of \f$ Z \f$ which represents the depth in the 3D camera frame.

  \param Z : \f$ Z \f$ value to set.
*/
void
CFeatureDefocusedLuminance::set_Z(const double Z)
{
    this->Z = Z ;
    flags[0] = true;
}


/*!
  Get the value of \f$ Z \f$ which represents the depth in the 3D camera frame.

  \return The value of \f$ Z \f$.
*/
double
CFeatureDefocusedLuminance::get_Z() const
{
    return Z ;
}


void
CFeatureDefocusedLuminance::setCameraParameters(CCameraThinLensParameters &_cam) 
{
  cam = _cam ;
}

void
CFeatureDefocusedLuminance::resetCameraParameters(CCameraThinLensParameters &_cam) 
{
  cam = _cam ;

  recompute_xy = true;
}



/*!

  Build a luminance feature directly from the image
*/
void
CFeatureDefocusedLuminance::buildFrom(vpImage<unsigned char> &I)
{
  double Ix,Iy, Ixx ;
  double px = cam.get_px(), py = cam.get_py();
	double x=0,y=0;
  CLuminanceDefocusCIP *pt_pixInfo;
  double *pt_s;
  int i, j;
  if (firstTimeIn==0)
  { 
    firstTimeIn=1 ;
    pt_pixInfo = (CLuminanceDefocusCIP *)pixInfo;
    pt_s = s.data;

    for(int n = 0 ; n < dim_s ; n++, pt_pixInfo++, pt_s++)
    {
      if(pt_pixInfo->toProcess)
      {
        pt_pixInfo->notToProcessCurrently = false;

        i = pt_pixInfo->i;
        j = pt_pixInfo->j;
        vpPixelMeterConversion::convertPoint(cam,
                                             j, i,
                                             x, y);


              Ix =  px * CImageFilter::derivativeFilterX(I, i, j);
              Iy =  py * CImageFilter::derivativeFilterY(I, i, j);
            
         Ixx = CImageFilter::laplacianFilterX(I, i, j) + CImageFilter::laplacianFilterY(I, i, j);

        pt_pixInfo->x = x;
        pt_pixInfo->y = y;
        pt_pixInfo->Z = Z;
        *pt_s  =  I[i][j] ;
        pt_pixInfo->I  = *pt_s;
        pt_pixInfo->Ix  = Ix;
        pt_pixInfo->Iy  = Iy;

        pt_pixInfo->Ixx  = Ixx;
      }
      else
        pt_pixInfo->notToProcessCurrently = true;
    }
  }
  else
  {
    pt_pixInfo = (CLuminanceDefocusCIP *)pixInfo;
    pt_s = s.data;
    for(int n = 0 ; n < dim_s ; n++, pt_pixInfo++, pt_s++)
    {
      if(pt_pixInfo->toProcess)
      {
        pt_pixInfo->notToProcessCurrently = false;

        i = pt_pixInfo->i;
        j = pt_pixInfo->j;


              Ix =  px * CImageFilter::derivativeFilterX(I, i, j);
              Iy =  py * CImageFilter::derivativeFilterY(I, i, j);

         Ixx = CImageFilter::laplacianFilterX(I, i, j) + CImageFilter::laplacianFilterY(I, i, j);

        *pt_s  =  I[i][j] ;
        pt_pixInfo->I  =  *pt_s;
        pt_pixInfo->Ix  = Ix;
        pt_pixInfo->Iy  = Iy;

        pt_pixInfo->Ixx  = Ixx;
      }
      else
        pt_pixInfo->notToProcessCurrently = true;
    }
  }
}



/*!

  Build a luminance feature directly from the image
*/

void
CFeatureDefocusedLuminance::buildFrom(vpImage<unsigned char> &I, CFeatureDefocusedLuminance *sd)
{
  if(sd == NULL)
    buildFrom(I);
  else
  {
   double Ix,Iy, Ixx ;
  double px = cam.get_px(), py = cam.get_py();
	double x=0,y=0;
  CLuminanceDefocusCIP *pt_pixInfo, *pt_pixInfod;
  double *pt_s;
  int i, j;//, cpt = 0;
  if (firstTimeIn==0)
  { 
    firstTimeIn=1 ;
    pt_pixInfo = (CLuminanceDefocusCIP *)pixInfo;
    pt_pixInfod = (CLuminanceDefocusCIP *)(sd->pixInfo);
    pt_s = s.data;

    for(int n = 0 ; n < dim_s ; n++, pt_pixInfo++, pt_pixInfod++, pt_s++)
    {
      pt_pixInfo->toProcess = pt_pixInfod->toProcess;
      pt_pixInfo->notToProcessCurrently = pt_pixInfod->notToProcessCurrently;
      if(pt_pixInfo->toProcess)
      {
	      i = pt_pixInfo->i;
	      j = pt_pixInfo->j;

        vpPixelMeterConversion::convertPoint(cam,
                                            j, i,
                                            x, y);

        pt_pixInfo->x = x;
        pt_pixInfo->y = y;
        pt_pixInfo->Z = Z;

        if(!pt_pixInfod->notToProcessCurrently)
        {
          
                Ix =  px * CImageFilter::derivativeFilterX(I, i, j);
                Iy =  py * CImageFilter::derivativeFilterY(I, i, j);
              
          Ixx = CImageFilter::laplacianFilterX(I, i, j) + CImageFilter::laplacianFilterY(I, i, j);

          *pt_s  =  I[i][j] ;
          pt_pixInfo->I  = *pt_s;
          pt_pixInfo->Ix  = Ix;
          pt_pixInfo->Iy  = Iy;

          pt_pixInfo->Ixx  = Ixx;
          //cpt++;
        }
      }
    }
  }
  else
  {
    pt_pixInfo = (CLuminanceDefocusCIP *)pixInfo;
    pt_pixInfod = (CLuminanceDefocusCIP *)(sd->pixInfo);
    pt_s = s.data;
    for(int n = 0 ; n < dim_s ; n++, pt_pixInfo++, pt_pixInfod++, pt_s++)
    {
      pt_pixInfo->toProcess = pt_pixInfod->toProcess;
      pt_pixInfo->notToProcessCurrently = pt_pixInfod->notToProcessCurrently;
      if(!pt_pixInfod->notToProcessCurrently)
      {
        i = pt_pixInfo->i;
        j = pt_pixInfo->j;

              Ix =  px * CImageFilter::derivativeFilterX(I, i, j);
              Iy =  py * CImageFilter::derivativeFilterY(I, i, j);

        Ixx = CImageFilter::laplacianFilterX(I, i, j) + CImageFilter::laplacianFilterY(I, i, j);

        *pt_s  =  I[i][j] ;
        pt_pixInfo->I  =  *pt_s;
        pt_pixInfo->Ix  = Ix;
        pt_pixInfo->Iy  = Iy;
        pt_pixInfo->Ixx  = Ixx;
        //cpt++;
      }
    }
  }
  //std::cout << "CFeatureLuminanceOmni::cartImagePlaneBuildFrom(I,sd) - cpt : " << cpt << std::endl;
  std::cout << "fin" << std::endl;

  }
}



/*!

  Compute and return the interaction matrix \f$ L_I \f$. The computation is made
  thanks to the values of the luminance features \f$ I \f$
*/
void
CFeatureDefocusedLuminance::interaction(vpMatrix &L)
{
  double *pt_L;
	CLuminanceDefocusCIP *pt_pixInfo = (CLuminanceDefocusCIP *)pixInfo;
  bool *pt_DOF;

  double f = cam.get_f();
  double x,y,Ix,Iy,Ixx;
  double xy, Z, iZ;
  int i;//, cpt = 0;

  L.resize(dim_s, nbDOF, false);
  pt_L = L.data;

  double cste = -cam.get_D()*f/(6.0*cam.get_ku()*(cam.get_Zf()-f)), IxxcsteilZ; //check -

  if(recompute_xy)
  {
    int i, j;
    recompute_xy = false;

	  for(int n = 0; n < dim_s; n++, pt_pixInfo++)
    {
		  if(pt_pixInfo->notToProcessCurrently)
		  {
			  for(i =0 ; i < nbDOF ; i++, pt_L++) *pt_L = 0. ;
		  }
		  else
		  {
        i = pt_pixInfo->i;
        j = pt_pixInfo->j;
        vpPixelMeterConversion::convertPoint(cam,
                                             j, i,
                                             x, y);
		    pt_pixInfo->x = x;
		    pt_pixInfo->y = y;

	      Ix = pt_pixInfo->Ix;
	      Iy = pt_pixInfo->Iy;
        Ixx = pt_pixInfo->Ixx;
        Z = pt_pixInfo->Z;

        iZ = 1./Z;
        IxxcsteilZ = Ixx*cste/(Z-f);

        xy = x*y;

        pt_DOF = DOF;
        if(*pt_DOF)     { *pt_L = -(Ix * (-iZ)); pt_L++; } // + Iy*0 + Ixx*cstesZ*0
        if(*(++pt_DOF)) { *pt_L = -(              Iy * (-iZ)); pt_L++; } // Ix*0 +  + Ixx*cstesZ*0
        if(*(++pt_DOF)) { *pt_L = -(Ix * (x*iZ) + Iy * (y*iZ)       + IxxcsteilZ*(-1)); pt_L++; }
        if(*(++pt_DOF)) { *pt_L = -(Ix * (xy)   + Iy * (1. + y*y)   + IxxcsteilZ*(-y*Z)); pt_L++; }
        if(*(++pt_DOF)) { *pt_L = -(Ix * (-(1. + x*x)) + Iy * (-xy) + IxxcsteilZ*x*Z); pt_L++; }
        if(*(++pt_DOF)) { *pt_L = -(Ix * (y)          + Iy * (-x)); pt_L++; } // + Ixx*cstesZ*0

          //cpt++;
		  }
    }
  }
  else
  {
	  for(int n = 0; n < dim_s; n++, pt_pixInfo++)
    {
		  if(pt_pixInfo->notToProcessCurrently)
		  {
			  for(i =0 ; i < nbDOF ; i++, pt_L++) *pt_L = 0. ;
		  }
		  else
		  {
		    x = pt_pixInfo->x;
		    y = pt_pixInfo->y;

	      Ix = pt_pixInfo->Ix;
	      Iy = pt_pixInfo->Iy;
        Ixx = pt_pixInfo->Ixx;
        Z = pt_pixInfo->Z;

        iZ = 1./Z;
        IxxcsteilZ = Ixx*cste/(Z-f);

        xy = x*y;

        pt_DOF = DOF;
        if(*pt_DOF)     { *pt_L = -(Ix * (-iZ)); pt_L++; } // + Iy*0 + Ixx*cstesZ*0
        if(*(++pt_DOF)) { *pt_L = -(              Iy * (-iZ)); pt_L++; } // Ix*0 +  + Ixx*cstesZ*0
        if(*(++pt_DOF)) { *pt_L = -(Ix * (x*iZ) + Iy * (y*iZ)       + IxxcsteilZ*(-1)); pt_L++; }
        if(*(++pt_DOF)) { *pt_L = -(Ix * (xy)   + Iy * (1. + y*y)   + IxxcsteilZ*(-y*Z)); pt_L++; }
        if(*(++pt_DOF)) { *pt_L = -(Ix * (-(1. + x*x)) + Iy * (-xy) + IxxcsteilZ*x*Z); pt_L++; }
        if(*(++pt_DOF)) { *pt_L = -(Ix * (y)          + Iy * (-x)); pt_L++; } // + Ixx*cstesZ*0

          //cpt++;
		  }
    }
  }
  //std::cout << "CFeatureLuminanceOmni::cartImagePlaneInteraction - cpt : " << cpt << std::endl;

}

/*!
  Compute and return the interaction matrix \f$ L_I \f$. The computation is made
  thanks to the values of the luminance features \f$ I \f$
*/
vpMatrix
CFeatureDefocusedLuminance::interaction(unsigned int /* select */)
{
  static vpMatrix L  ;
  interaction(L);
  return L ;
}


/*!
  Compute the error \f$ (I-I^*)\f$ between the current and the desired
 
  \param s_star : Desired visual feature.
  \param e : Error between the current and the desired features.

*/
void
CFeatureDefocusedLuminance::error(const vpBasicFeature &s_star,
			  vpColVector &e)
{
  e.resize(dim_s) ;
  double *pt_e = e.data, *pt_s = s.data;
  //int cpt = 0;
  
        CLuminanceDefocusCIP *pt_pixInfo = (CLuminanceDefocusCIP *)pixInfo;
        for (int i =0 ; i < dim_s ; i++, pt_pixInfo++, pt_e++, pt_s++)
        {
          if(pt_pixInfo->notToProcessCurrently)
            *pt_e = 0.;
          else
          {
            *pt_e = *pt_s - s_star[i] ;
            //cpt++;
          }
        }
  //std::cout << "CFeatureLuminanceOmni::error - cpt : " << cpt << std::endl;
}



/*!
  Compute the error \f$ (I-I^*)\f$ between the current and the desired
 
  \param s_star : Desired visual feature.
  \param select : Not used.

*/
vpColVector
CFeatureDefocusedLuminance::error(const vpBasicFeature &s_star,
			  const int /* select */)
{
  static vpColVector e ;
  
  error(s_star, e) ;
  
  return e ;

}

/*!

  Not implemented.

 */
void
CFeatureDefocusedLuminance::print(unsigned int /* select */) const
{
  static int firsttime =0 ;

  if (firsttime==0)
  {
    firsttime=1 ;
    vpERROR_TRACE("not implemented") ;
    // Do not throw and error since it is not subject
    // to produce a failure
  }
 }



/*!

  Not implemented.

 */
void
CFeatureDefocusedLuminance::display(const vpCameraParameters & /* cam */,
			    const vpImage<unsigned char> & /* I */,
			    const vpColor & /* color */,  unsigned int /* thickness */) const
{
 static int firsttime =0 ;

  if (firsttime==0)
  {
    firsttime=1 ;
    vpERROR_TRACE("not implemented") ;
    // Do not throw and error since it is not subject
    // to produce a failure
  }
}

/*!

  Not implemented.

 */
void
CFeatureDefocusedLuminance::display(const vpCameraParameters & /* cam */,
			    const vpImage<vpRGBa> & /* I */,
			    const vpColor & /* color */, unsigned int /* thickness */) const
{
  static int firsttime =0 ;

  if (firsttime==0)
  {
    firsttime=1 ;
    vpERROR_TRACE("not implemented") ;
    // Do not throw and error since it is not subject
    // to produce a failure
  }
}


/*!
  Create an object with the same type.

  \code
  vpBasicFeature *s_star;
  vpFeatureLuminance s;
  s_star = s.duplicate(); // s_star is now a vpFeatureLuminance
  \endcode

*/
CFeatureDefocusedLuminance*
CFeatureDefocusedLuminance::duplicate() const
{
  CFeatureDefocusedLuminance *feature = new CFeatureDefocusedLuminance ;
  return feature ;
}

/*
 * Local variables:
 * c-basic-offset: 2
 * End:
 */
