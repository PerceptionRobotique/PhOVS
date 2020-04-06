/****************************************************************************
 *
 * July 2011
 *
 * Author:
 * Guillaume Caron
 * enriching from visp/vpImageTools
 *
 *****************************************************************************/

#ifndef CImageTools_H
#define CImageTools_H

/*!
  \file CImageTools.h
  \brief  Various image tools (ZNCC, ...)

*/

#include <fstream>
#include <iostream>
#include <math.h>
#include <string.h>

#include <visp/vpConfig.h>
#include <visp/vpImageException.h>
#include <visp/vpImage.h>
#include <visp/vpMatrix.h>
#include <visp/vpMath.h>
#include <visp/vpImageTools.h>

/*!
  \class CImageTools

  \ingroup ImageTool

  \brief  Various image tools (ZNCC, ...)

*/
class CImageTools : public vpImageTools
{

public:
  static void imageDifference(vpImage<unsigned char> &I1, 
			                        vpImage<unsigned char> &I2,
			                        vpImage<unsigned char> &Idiff) ; //surcharge : parcours efficace des images

  static void imageDifference(vpImage<unsigned char> &I1, 
				   vpImage<unsigned char> &I2,
           vpImage<float> &IDepthMask,
				   vpImage<unsigned char> &Idiff);

  static void imageAugment(vpImage<unsigned char> &I1, 
				   vpImage<unsigned char> &I2,
           vpImage<float> &IDepthMask,
				   vpImage<unsigned char> &IAugmented,
           float alpha = 1.0f);

  static double imageZNCC(vpImage<unsigned char> &I1, 
			                    vpImage<unsigned char> &I2) ;

  static void imageZN(vpImage<unsigned char> &I1, 
			                vpColVector &V);

  static void imageZN(vpImage<unsigned char> &I1, 
			                vpColVector &V,
                      vpImage<float> &Idepth_Mask,
                          bool ZN_masked_only = false);

  static void imageZN(vpImage<unsigned char> &I1, 
			                vpImage<double> &I2);

  static void imageZN(vpImage<unsigned char> &I1, 
			                vpImage<double> &I2,
                      vpImage<float> &Idepth_Mask,
                          bool ZN_masked_only = false);

  static double imageZNCC(vpImage<unsigned char> &I1, 
			                    vpImage<double> &I2);

  static double imageZNCC(vpImage<unsigned char> &I1, 
			                    vpColVector &vI2zn);

} ;


#endif


/*
 * Local variables:
 * c-basic-offset: 2
 * End:
 */
