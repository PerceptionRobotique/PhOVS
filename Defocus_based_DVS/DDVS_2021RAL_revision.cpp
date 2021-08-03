//#define TESTCOMPCROBOT
//#define TESTACQUICAM

//./DDVS 1000 0 400 -400 27 0 0
// factor to meters
// tX
// tY
// tZ
// rX
// rY
// rZ

#include "src/C_UR.h"

#include <PhOVS/CImageTools.h>
#include <PhOVS/CCameraThinLensParameters.h>
//#include <PhOVS/CFeatureLuminanceOmni.h>
#include <PhOVS/CFeatureDefocusedLuminance.h>

#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpDisplayX.h>

#ifndef TESTCOMPCROBOT
	#include "src/CamFlir.hpp"
#endif


//Pour images acquises de 640x512
#define ACQWIDTH 640
#define ACQHEIGHT 512
#define DI 5
#define DJ 5
#define FACT 0.25

/*
//Pour images acquises de 1024x1024
#define ACQWIDTH 1024
#define ACQHEIGHT 1024
#define DI 5
#define DJ 5
#define FACT 0.25
*/

/*
//Pour images acquises de 640x492 (mais 616x492 desires)
#define ACQWIDTH 640
#define ACQHEIGHT 492
#define DI 5
#define DJ (12+5)
#define FACT 0.5
*/

/*
//Pour images acquises de 704x684 (mais 684x484 desires)
#define ACQWIDTH 704
#define ACQHEIGHT 684
#define DI 5
#define DJ (10+5)
#define FACT 0.5
*/


#define NBRI (ACQHEIGHT*FACT-2*DI)
#define NBRJ (ACQWIDTH*FACT-2*DJ)

//#define PAS (4*FACT)

#define PAS 1

#define Z 0.75
//0.5


#define CURRENT

//#define FACT_DDL_NON_VZ 1.5
#define FACT_DDL_NON_VZ 1.0

#define INDICATORS
#define FILE_EXT "jpg"
#define OPT_DISP
//#define OPT_CLICK


int
main(int argc, const char ** argv)
{
#ifdef INDICATORS
    std::ostringstream s;
    std::string filename;
#endif // INDICATORS

#ifndef TESTACQUICAM
	C_UR UR10("192.168.1.3", 30003, 0.12); //30002: 10 Hz ; 30003: 125/500 Hz

  // with F0.95 lens
	vpColVector j_init(6);
  //Flir FL3
/*
  //desk, Marylin
	j_init[0] = vpMath::rad(-51.13);
	j_init[1] = vpMath::rad(-129.05);
	j_init[2] = vpMath::rad(-96.63);
	j_init[3] = vpMath::rad(-41.79);
	j_init[4] = vpMath::rad(88.69);
	j_init[5] = vpMath::rad(84.95);
*/
/*
  //red table
  //0.25 m depth
	j_init[0] = vpMath::rad(-118.10);
	j_init[1] = vpMath::rad(-100.44);
	j_init[2] = vpMath::rad(-137.12);
	j_init[3] = vpMath::rad(-30.09);
	j_init[4] = vpMath::rad(88.61);
	j_init[5] = vpMath::rad(94.94);
*/
/*
  //objects to grasp
  //0.25 m depth
	j_init[0] = vpMath::rad(-170.17);
	j_init[1] = vpMath::rad(-159.73);
	j_init[2] = vpMath::rad(-112.02);
	j_init[3] = vpMath::rad(94.55);
	j_init[4] = vpMath::rad(85.74);
	j_init[5] = vpMath::rad(88.81);
*/
/*
  //ground
  //0.25 m depth
	j_init[0] = vpMath::rad(-135.44);
	j_init[1] = vpMath::rad(-186.68);
	j_init[2] = vpMath::rad(-72.01);
	j_init[3] = vpMath::rad(-6.46);
	j_init[4] = vpMath::rad(90.63);
	j_init[5] = vpMath::rad(91.66);
*/
/*
  //ground
  //0.5 m depth
	j_init[0] = vpMath::rad(-135.91);
	j_init[1] = vpMath::rad(-163.16);
	j_init[2] = vpMath::rad(-95.35);
	j_init[3] = vpMath::rad(-6.66);
	j_init[4] = vpMath::rad(90.63);
	j_init[5] = vpMath::rad(91.46);
*/
/*
  //ground
  //0.5 m depth 2
	j_init[0] = vpMath::rad(-130.02);
	j_init[1] = vpMath::rad(-159.33);
	j_init[2] = vpMath::rad(-99.84);
	j_init[3] = vpMath::rad(-5.89);
	j_init[4] = vpMath::rad(90.14);
	j_init[5] = vpMath::rad(97.37);
*/

/*
  //Flir GS3
  //ground
  //0.25 m depth
	j_init[0] = vpMath::rad(-135.47);
	j_init[1] = vpMath::rad(-184.98);
	j_init[2] = vpMath::rad(-74.05);
	j_init[3] = vpMath::rad(-6.12);
	j_init[4] = vpMath::rad(90.63);
	j_init[5] = vpMath::rad(91.64);
*/

  //ground
  //0.5 m depth 2
	j_init[0] = vpMath::rad(-130.02);
	j_init[1] = vpMath::rad(-159.33);
	j_init[2] = vpMath::rad(-99.84);
	j_init[3] = vpMath::rad(-5.89);
	j_init[4] = vpMath::rad(90.14);
	j_init[5] = vpMath::rad(97.37);

  UR10.setCameraArticularPose(j_init);

  vpTime::wait(3000);

/*	vpColVector v(6);
	v[0] = -0.1;
	v[1] = 0;
	v[2] = 0;
	v[3] = 0;
	v[4] = 0;
	v[5] = 0.0;
  UR10.setCameraVelocity(v);
  vpTime::wait(2000);
	v[0] = 0;
	v[1] = 0;
	v[2] = 0;
	v[3] = 0;
	v[4] = 0;
	v[5] = -0.5;
  UR10.setCameraVelocity(v);
  vpTime::wait(2000);
	v[0] = 0.1;
	v[1] = 0;
	v[2] = 0;
	v[3] = 0;
	v[4] = 0;
	v[5] = 0.0;
  UR10.setCameraVelocity(v);
  vpTime::wait(2000);

	v[0] = 0;
	v[1] = 0;
	v[2] = 0;
	v[3] = 0.5;
	v[4] = 0;
	v[5] = 0.0;
  UR10.setCameraVelocity(v);
  vpTime::wait(2000);
	v[0] = 0;
	v[1] = 0;
	v[2] = -0.5;
	v[3] = 0;
	v[4] = 0;
	v[5] = 0.0;
  UR10.setCameraVelocity(v);
  vpTime::wait(2000);
  v[0] = 0;
	v[1] = 0;
	v[2] = 0;
	v[3] = 0.0;
	v[4] = 0.5;
	v[5] = 0.0;
  UR10.setCameraVelocity(v);
  vpTime::wait(2000);
	v[0] = 0;
	v[1] = 0;
	v[2] = 0.5;
	v[3] = 0;
	v[4] = 0.0;
	v[5] = 0.0;
  UR10.setCameraVelocity(v);
  vpTime::wait(2000);
*/

#endif
#ifndef TESTCOMPCROBOT
	vpImage<unsigned char> Iacq, I, Id;

  
  //Parametres intrinseques pour FlirCam FL3 1/1.8" 1280x1024 5.3 um (matrice: 6.784 mm x 5.4272 mm)
//	int larg = ACQWIDTH*FACT, haut = ACQHEIGHT*FACT;
  // Fujinon lens
//  double px = 500*FACT;//1200;
//  double py = 500*FACT;//1200;
  // F0.95 lens
  //double px = 1603.725*FACT;//1200;
  //double py = 1603.725*FACT;//1200;
  double ku = (2*5.3e-6)/FACT; // m //ATTENTION: considere pour resultats RAL-ICRA 2021 a la place de (10.6e-6)/FACT
  

  //Parametres intrinseques pour FlirCam GS3 1" 2048x2048 5.5 um (matrice: 11.264 mm x 11.264 mm)
  
  //Maxi pour le capteur
	int larg = ACQWIDTH*FACT, haut = ACQHEIGHT*FACT;
  
  /*
  //Simuler la FL3
	//int larg = (1024*6.784/11.264)*FACT, haut = (1024*5.4272/11.264)*FACT;
  int larg = ACQWIDTH*FACT, haut = ACQHEIGHT*FACT;
  */  
  /*
  //Maxi pour l'objectif Yakumo F0.95 (2/3")
	int larg = ACQWIDTH*FACT, haut = ACQHEIGHT*FACT;
  */
  double u0 = (ACQWIDTH*0.5)*FACT;
  double v0 = (ACQHEIGHT*0.5)*FACT;  
  
  // Fujinon lens
//  double px = 500*FACT;//1200;
//  double py = 500*FACT;//1200;
  // F0.95 lens
  //double px = 1603.725*FACT;//1200;
  //double py = 1603.725*FACT;//1200;  
  //double ku = (11.0e-6)/FACT; // m


  //Parametres objectif Yakumo
  double precond = 0.1;//1.0;//
  //ku /= precond;
  double f = precond*17e-3; // m
  double FNumber = 0.95; //no unit
  double Zf = 0.5;//0.25; // m

  CCameraThinLensParameters ccam(f, ku, FNumber, Zf, u0, v0);

  CamFlir<unsigned char> grabber(ACQWIDTH,ACQHEIGHT,8,0); 
 
  //Acquisition
  grabber.getFrame(Iacq);

  vpImageTools::resize(Iacq, I, larg, haut, vpImageTools::vpImageInterpolationType::INTERPOLATION_CUBIC);//INTERPOLATION_LINEAR);//

  std::cout << I.getHeight() << " " << I.getWidth() << std::endl;

  Id = I;

  // display the desired image
  #if defined VISP_HAVE_X11
  vpDisplayX dId;
  #elif defined VISP_HAVE_GDI
  vpDisplayGDI dId;
  #elif defined VISP_HAVE_GTK
  vpDisplayGTK dId;
  #endif

#if defined(VISP_HAVE_X11) || defined(VISP_HAVE_GDI) || defined(VISP_HAVE_GTK) 
//  if (opt_display) 
  {
    dId.init(Id, 30, 10, "Photometric central visual servoing : s*") ;
  }
#endif

  vpDisplay::display ( Id ) ;
  vpDisplay::flush ( Id ) ;
#ifdef OPT_CLICK
  vpDisplay::getClick(Id);
#endif //OPT_CLICK

#ifdef INDICATORS
  s.str("");
  s.setf(std::ios::right, std::ios::adjustfield);
  s << "resultat/Id." << FILE_EXT;
  filename = s.str();
  vpImageIo::write(Id, filename);

  //save the desired pose to file
  s.str("");
  s.setf(std::ios::right, std::ios::adjustfield);
  s << "resultat/desiredPose.txt";
  filename = s.str();
  std::ofstream ficDesiredPose(filename.c_str());

  vpColVector p;
  UR10.getCameraPoseRaw(p);
  ficDesiredPose << p.t() << std::endl;

  ficDesiredPose.close();
#endif //INDICATORS
  // ------------------------------------------------------
  // Visual feature, interaction matrix, error
  // s, Ls, Lsd, Lt, Lp, etc
  // ------------------------------------------------------

  // desired visual feature built from the image
  CFeatureDefocusedLuminance sId ;
  sId.setCameraParameters(ccam) ;
  sId.init(Id.getHeight(), Id.getWidth(), DI, DJ, NBRI, NBRJ, PAS, NULL, Z) ;
  bool dofs[6] = {true, true, true, true, true, true}; // 6 DoF
//  bool dofs[6] = {true, true, true, false, false, true}; // 4 DoF
//  bool dofs[6] = {false, true, true, true, false, false};
  sId.set_DOF(dofs[0], dofs[1], dofs[2], dofs[3], dofs[4], dofs[5]);
std::cout << "ok 1" << std::endl;
  sId.buildFrom(Id);
std::cout << "ok 2" << std::endl;
  // Matrice d'interaction, Hessien, erreur,...
  vpMatrix Lsd, piLsd;   // matrice d'interaction a la position desiree
  vpColVector error ; // Erreur I-I*

#ifndef CURRENT
  // Compute the interaction matrix
  // link the variation of image intensity to camera motion
  // here it is computed at the desired position
  sId.interaction(Lsd) ;
  piLsd = Lsd.pseudoInverseEigen3();
#endif
std::cout << "ok 3" << std::endl;


  std::cout << "Deplacement vers pose initiale " << argc << " " << atof(argv[2])/atof(argv[1]) << " " << atof(argv[3])/atof(argv[1]) << std::endl;
	vpColVector p_init;
  p_init.resize(6);

  //argv[1] : metric scale
  //argv[2] : tX
  //argv[3] : tY
  //argv[4] : tZ
  //argv[5] : rX directly expected in degrees
  //argv[6] : rY directly expected in degrees
  //argv[7] : rZ directly expected in degrees

  if(argc > 2)
  {
    p_init[0] = atof(argv[2])/atof(argv[1]);

    std::cout << p_init.t() << std::endl;

    if(argc > 3)
    {
      p_init[1] = atoi(argv[3])/atof(argv[1]);

      if(argc > 4)
      {
        p_init[2] = atoi(argv[4])/atof(argv[1]);

        if(argc > 5)
        {
          p_init[3] = atoi(argv[5])*M_PI/180.0;

          if(argc > 6)
          {
            p_init[4] = atoi(argv[6])*M_PI/180.0;

            if(argc > 7)
            {
              p_init[5] = atoi(argv[7])*M_PI/180.0;
            }
          }
        }
      }
    }


    std::cout << p_init.t() << std::endl;
  }
/*	p_init[0] = 0.0;
	p_init[1] = 0;
	p_init[2] = -0.25; //-0.06
	p_init[3] = 0;
	p_init[4] = 0;
	p_init[5] = 0.0;*/

	UR10.setCameraRelativePose(p_init);

  vpTime::wait(7000);

  //Acquisition
  grabber.getFrame(Iacq);
  vpImageTools::resize(Iacq, I, larg, haut, vpImageTools::vpImageInterpolationType::INTERPOLATION_CUBIC);//INTERPOLATION_LINEAR);//

   // ------------------------------------------------------
    // Visual feature, interaction matrix, error
    // s, Ls, Lsd, Lt, Lp, etc
    // ------------------------------------------------------
   
   // current visual feature built from the image 
    // (actually, this is the image...)
    CFeatureDefocusedLuminance sI ;
    sI.setCameraParameters(ccam) ;
    sI.init(I.getHeight(), I.getWidth(), DI, DJ, NBRI, NBRJ, PAS, NULL, Z) ;
    sI.set_DOF(dofs[0], dofs[1], dofs[2], dofs[3], dofs[4], dofs[5]);

    sI.buildFrom(I, &sId);

	// initialise a  display
  vpDisplayX display;

  display.init(I, 100, 100, "Photometric central visual servoing : s");

  vpDisplay::display(I);
  vpDisplay::flush(I);

  vpImage<unsigned char> Idiff(haut, larg) ;

  vpImageTools::imageDifference(I,Id,Idiff) ;

  // Affiche de l'image de difference
  #if defined VISP_HAVE_X11
  vpDisplayX d1;
  #elif defined VISP_HAVE_GDI
  vpDisplayGDI d1;
  #elif defined VISP_HAVE_GTK
  vpDisplayGTK d1;
  #endif

#if defined(VISP_HAVE_X11) || defined(VISP_HAVE_GDI) || defined(VISP_HAVE_GTK) 
//  if (opt_display) {
    d1.init(Idiff, 680, 10, "Photometric central visual servoing : s-s* ") ;
//  }
#endif

    vpDisplay::display(Idiff);
    vpDisplay::flush(Idiff);
#ifdef OPT_CLICK
    vpDisplay::getClick(Idiff);
#endif // OPT_CLICK
    // ------------------------------------------------------
    // Control law
    double lambda ; //gain
    vpColVector e ;
    vpColVector v ; // camera velocity sent to the robot

    //lambda   = 4.0;//2.8;//Fmax 4 DoF
    //lambda   = 1.0;//2.8;//Fmax 6 DoF

    //lambda   = 4.0;//2.8;//Fmin 4 DoF //Marylin
    //lambda   = 3.0;//2.8;//Fmin 4 DoF //objects
    //lambda   = 2.0;//Fmin ; // 6 ddl precond, Z 50 cm
    //lambda   = 0.5;//0.5;//Fmin ; // 6 ddl precond, Z 25 cm
    //lambda   = 0.5;//Fmin ; // 6 ddl sans precond

    //lambda   = 4.0;//2.8;//Fmin 4 DoF //Marylin
    lambda   = 1.0;//0.5;//Fmin ; // 6 ddl //Marylin

    // ----------------------------------------------------------
    int iter   = 1;
vpMatrix Ls, piLs;

  unsigned int nbDOF = 6, numDOF, indDOF;
  vpColVector v6(6);
  double residual;
#ifdef INDICATORS
/*  vpPoseVector p;
  std::vector<vpPoseVector> v_p;*/
  std::vector<vpColVector> v_p;
  std::vector<double> v_residuals;
  std::vector<vpImage<unsigned char> > v_I;
  std::vector<vpImage<unsigned char> > v_Idiff;
  std::vector<double> v_tms;
  std::vector<bool> v_servo_actif;
  vpMatrix V;
  vpColVector w;
  std::vector<double> v_svd;
  double cond, cond_back;
#endif //INDICATORS
  double tms, duree;
  vpMouseButton::vpMouseButtonType btn;
  bool btn_pressed = false;
  bool servo_actif = false;

/*      ccam.set_f(f * 0.1);
      sI.resetCameraParameters(ccam) ;
*/
  bool first_part = true;
	do
	{
      std::cout << "--------------------------------------------" << iter++ << std::endl ;

#ifdef INDICATORS
    UR10.getCameraPoseRaw(p);
#endif //INDICATORS

		grabber.getFrame(Iacq);

      tms = vpTime::measureTimeMs();
    vpImageTools::resize(Iacq, I, larg, haut, vpImageTools::vpImageInterpolationType::INTERPOLATION_CUBIC);//INTERPOLATION_LINEAR);//

    sI.buildFrom(I, &sId) ;

      // compute current error
      sI.error(sId,error) ;

      residual = 0.5*error.sumSquare();
      std::cout << "error : " << residual << std::endl;

      /*
      if(first_part && (residual < 5e7))
      {
        first_part = false;
        ccam.set_f(f);
        sI.resetCameraParameters(ccam) ;
      }*/

      if(residual < 1e7)
        lambda=0.5;

/*      if(iter > 60)
      {
              ccam.set_f(f);
          sI.resetCameraParameters(ccam) ;
      }
*/
#ifdef CURRENT
/*
    if( (iter > 2) )
    {
      precond *= 10.0/cond;
      if(precond > 10) precond = 10;
      if(precond < 0.1) precond = 0.1;

      ccam.set_f(f * precond);
      sI.resetCameraParameters(ccam) ;
    }
*/
    sI.interaction(Ls) ;
    piLs = Ls.pseudoInverseEigen3();

      //  compute the control law 
      e = piLs*error;
#else
      e = piLsd*error;
#endif
      v = - lambda*e;

      //v6 = v;//0;
      //v6[0] = v[0]*FACT_DDL_NON_VZ;
      //v6[2] = v[1];
      //v6[4] = v[0];
/*
    std::cout << "v : " << v.t() << std::endl;

    double normv = sqrt(v.sumSquare());

    std::cout << "normv : " << normv << std::endl;

    if(normv > 0.1*lambda)
      v *= 0.1*lambda/normv;
*/
    //update the DOFs
    indDOF = 0;
    for (numDOF = 0 ; numDOF < nbDOF ; numDOF++)
        if (dofs[numDOF])
        {
            v6[numDOF] = v[indDOF];
            indDOF++;
        }
        else
            v6[numDOF] = 0;

      std::cout << "v6 : " << v6.t() << std::endl;

      if(servo_actif)
        UR10.setCameraVelocity(v6);

    duree = vpTime::measureTimeMs() - tms;
    std::cout << "duration : " << duree <<std::endl;

#if defined(OPT_DISP) || defined(INDICATORS)
      vpImageTools::imageDifference(I,Id,Idiff) ;
#endif

#ifdef OPT_DISP
		vpDisplay::display(I);
		vpDisplay::flush(I);

    vpDisplay::display(Idiff);
    vpDisplay::flush(Idiff);
#endif

#ifdef INDICATORS
    v_p.push_back(p);
    v_residuals.push_back(residual);
    v_I.push_back(I);
    v_Idiff.push_back(Idiff);
    v_tms.push_back(duree);
    v_servo_actif.push_back(servo_actif);

    Ls.svdEigen3(w, V);
    cond_back = cond;
    cond = w[0]/w[w.getRows()-1];
    v_svd.push_back(cond);
#endif //INDICATORS      
  
    btn_pressed=vpDisplay::getClick(I, btn, false);

    if(btn_pressed && (btn == vpMouseButton::button3))
      servo_actif = !servo_actif;

	}
	while(!btn_pressed || (btn != vpMouseButton::button1));

  v6.resize(6);
  UR10.setCameraVelocity(v6);

#ifdef INDICATORS
    //save pose list to file
    s.str("");
    s.setf(std::ios::right, std::ios::adjustfield);
    s << "resultat/poses.txt";
    filename = s.str();
    std::ofstream ficPoses(filename.c_str());
    //save residual list to file
    s.str("");
    s.setf(std::ios::right, std::ios::adjustfield);
    s << "resultat/residuals.txt";
    filename = s.str();
    std::ofstream ficResiduals(filename.c_str());
    //save the processing times
    s.str("");
    s.setf(std::ios::right, std::ios::adjustfield);
    s << "resultat/times.txt";
    filename = s.str();
    std::ofstream ficTimes(filename.c_str());
    //save svd list to file
    s.str("");
    s.setf(std::ios::right, std::ios::adjustfield);
    s << "resultat/cond.txt";
    filename = s.str();
    std::ofstream ficSVD(filename.c_str());
    //save servo actif list to file
    s.str("");
    s.setf(std::ios::right, std::ios::adjustfield);
    s << "resultat/servoActif.txt";
    filename = s.str();
    std::ofstream ficServo(filename.c_str());

    for(unsigned int i = 0 ; i < v_p.size() ; i++)
    {
      ficPoses << v_p[i].t() << std::endl;
      ficResiduals << std::fixed << std::setw( 11 ) << std::setprecision( 6 ) << v_residuals[i] << std::endl;
      ficTimes << v_tms[i] << std::endl;
      ficSVD << v_svd[i] << std::endl;
      ficServo << v_servo_actif[i] << std::endl;

      s.str("");
      s.setf(std::ios::right, std::ios::adjustfield);
      s << "resultat/I/I." << std::setw(4) << std::setfill('0') << i << "." << FILE_EXT;
      filename = s.str();
      vpImageIo::write(v_I[i], filename);

      s.str("");
      s.setf(std::ios::right, std::ios::adjustfield);
      s << "resultat/IErr/IErr." << std::setw(4) << std::setfill('0') << i << "." << FILE_EXT;
      filename = s.str();
      vpImageIo::write(v_Idiff[i], filename);
    }

    ficPoses.close();
    ficResiduals.close();
    ficTimes.close();
    ficSVD.close();
    ficServo.close();


#endif //INDICATORS      

#endif

	return 0;
}
