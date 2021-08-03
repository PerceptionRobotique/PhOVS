/****************************************************************************
 *
 * October 2020
 *
 * Author:
 * Guillaume Caron
 * inspired from CFeatureLuminanceOmni
 *
 *****************************************************************************/

#ifndef CFeatureDefocusedLuminance_h
#define CFeatureDefocusedLuminance_h

#include <visp/vpConfig.h>
#include <visp/vpMatrix.h>
#include <visp/vpBasicFeature.h>
#include <visp/vpImage.h>

#include <PhOVS/CCameraThinLensParameters.h>

/*!
  \class CLuminanceDefocus
  \brief Class that defines the brightness, gradient and laplacian of a point

  \sa CLuminanceDefocus
*/


class CLuminanceDefocus
{
 public:
  int i, j;
  double I ; // pixel intensity
  double Z; // pixel depth
  bool toProcess, notToProcessCurrently; // boolean to precise if this pixel is to be processed or not (low gradient, M-Estimator, etc)

  CLuminanceDefocus() { toProcess = true; notToProcessCurrently = false;}
};

//luminance in cartesian image plane
class CLuminanceDefocusCIP : public CLuminanceDefocus
{
 public:
  double x, y;   // point coordinates (in meter)
  double Ix,Iy ; // pixel gradient
  double Ixx ; // pixel laplacian
};


/*!
  \class CFeatureDefocusedLuminance
  \brief Class that defines the defocused image brightness visual feature

  for more details see
  G. Caron, RAL/ICRA 2021 submission
*/

class CFeatureDefocusedLuminance : public vpBasicFeature
{
public:

 protected:
  //! FeaturePoint depth (required to compute the interaction matrix)
  //! default Z = 1m
  double Z ;
  bool recompute_xy;

  int imWidth, imHeight, di, nbri, dj, nbrj, pas;
  int derivativeMaskHalfSize, nbNeigh;
  int nbDim;

  int nbDOF;
  bool DOF[6];
  
  //! Store the image (as a vector with intensity and gradient I, Ix, Iy) 
  CLuminanceDefocus *pixInfo ;
  int  firstTimeIn  ;

 public:
  void buildFrom(vpImage<unsigned char> &I) ;
  void buildFrom(vpImage<unsigned char> &I, CFeatureDefocusedLuminance *sd) ;

public: 
  void init() ;
  void init(int _imHeight, int _imWidth, int _di = 10, int _dj = 10, int _nbri = 10, int _nbrj = 10, int _pas = 1, vpImage<unsigned char> *Imask = NULL, double _rho = 1., int _derivativeMaskSize = 7) ;

  CFeatureDefocusedLuminance() ;
 
  //! Destructor.
  virtual ~CFeatureDefocusedLuminance() ;

 public:
  CCameraThinLensParameters cam ;
  void setCameraParameters(CCameraThinLensParameters &_cam) ;
  void resetCameraParameters(CCameraThinLensParameters &_cam) ;
  /*
    section Set/get Z
  */

  void set_DOF(bool un = true, bool deux = true, bool trois = true, bool quatre = true, bool cinq = true, bool six = true);
  void set_Z(const double Z) ;
  double get_Z() const  ;


  /*
    vpBasicFeature method instantiation
  */

 
  vpMatrix  interaction(unsigned int select = FEATURE_ALL);
  void      interaction(vpMatrix &L);
private:
  void cartImagePlaneInteraction(vpMatrix &L);
public:
  vpColVector error(const vpBasicFeature &s_star,
		    const int select = FEATURE_ALL)  ;
  void error(const vpBasicFeature &s_star,
	     vpColVector &e)  ;

  void print(unsigned int select = FEATURE_ALL ) const ;

  CFeatureDefocusedLuminance *duplicate() const ;


  void display(const vpCameraParameters &cam,
	       const vpImage<unsigned char> &I,
	       const vpColor &color, unsigned int thickness=1) const ;
  void display(const vpCameraParameters &cam,
               const vpImage<vpRGBa> &I,
               const vpColor &color, unsigned int thickness=1) const ;

  //! Compute the error between a visual features and zero
  vpColVector error(const int select = FEATURE_ALL)  ;

} ;

#endif //CFeatureDefocusedLuminance_h
