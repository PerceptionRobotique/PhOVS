/****************************************************************************
 *
 * July 2011
 *
 * Author:
 * Guillaume Caron
 * inspired from visp/vpImageFilter (Eric Marchand)
 *
 *****************************************************************************/

#include <PhOVS/CImageTools.h>

/*!
  Compute the signed difference between the two images I1 and I2 for 
  visualization issue : Idiff = I1-I2

  - pixels with a null difference are set to 128. 
  - A negative difference implies a pixel value < 128
  - A positive difference implies a pixel value > 128
  
  \param I1 : The first image.
  \param I2 : The second image.
  \param Idiff : The result of the difference.
*/
void CImageTools::imageDifference(vpImage<unsigned char> &I1, 
				   vpImage<unsigned char> &I2,
				   vpImage<unsigned char> &Idiff)
{
  if ((I1.getHeight() != I2.getHeight()) || (I1.getWidth() != I2.getWidth()))
  {
    throw (vpException(vpException::dimensionError, "The two images have not the same size"));
  }
  
  if ((I1.getHeight() != Idiff.getHeight()) || (I1.getWidth() != Idiff.getWidth()))
    Idiff.resize(I1.getHeight(), I1.getWidth(), false);
  
  int n = I1.getHeight() * I1.getWidth() ;
  int diff ;
  unsigned char *pt_I1 = I1.bitmap, *pt_I2 = I2.bitmap, *pt_Idiff = Idiff.bitmap;
  for (int b = 0; b < n ; b++, pt_I1++, pt_I2++, pt_Idiff++)
    {
      diff = *pt_I1 - *pt_I2 + 128;
      *pt_Idiff = (unsigned char) (vpMath::maximum(vpMath::minimum(diff, 255), 0));
    }
}

/*!
  Compute the signed difference between the two images I1 and I2, only 
  for non empty pixels for visualization issue : Idiff = I1-I2

  - pixels with a null difference are set to 128. 
  - A negative difference implies a pixel value < 128
  - A positive difference implies a pixel value > 128
  
  \param I1 : The first image.
  \param I2 : The second image.
  \param IDepthMask : The depth mask image.
  \param Idiff : The result of the difference.
*/
void CImageTools::imageDifference(vpImage<unsigned char> &I1, 
				   vpImage<unsigned char> &I2,
           vpImage<float> &IDepthMask,
				   vpImage<unsigned char> &Idiff)
{
  if ((I1.getHeight() != I2.getHeight()) || (I1.getWidth() != I2.getWidth()) || (I1.getHeight() != IDepthMask.getHeight()) || (I1.getWidth() != IDepthMask.getWidth()))
  {
    throw (vpException(vpException::dimensionError, "Images have not the same size"));
  }
  
  if ((I1.getHeight() != Idiff.getHeight()) || (I1.getWidth() != Idiff.getWidth()))
    Idiff.resize(I1.getHeight(), I1.getWidth(), false);
  
  int n = I1.getHeight() * I1.getWidth() ;
  int diff ;
  unsigned char *pt_I1 = I1.bitmap, *pt_I2 = I2.bitmap, *pt_Idiff = Idiff.bitmap;
  float *pt_IDepthMask = IDepthMask.bitmap;
  for (int b = 0; b < n ; b++, pt_I1++, pt_I2++, pt_Idiff++, pt_IDepthMask++)
    {
      if(*pt_IDepthMask != -1.f)
      {
        diff = *pt_I1 - *pt_I2 + 128;
        *pt_Idiff = (unsigned char) (vpMath::maximum(vpMath::minimum(diff, 255), 0));
      }
      else
        *pt_Idiff = *pt_I2;
    }
}

/*!
  Compute the blend of the two images I1 and I2 (background), only 
  for non empty pixels
  
  \param I1 : The first image.
  \param I2 : The second background image.
  \param IDepthMask : The depth mask image.
  \param Idiff : The result of the difference.
  \alpha the alpha blending coefficient
*/
void CImageTools::imageAugment(vpImage<unsigned char> &I1, 
				   vpImage<unsigned char> &I2,
           vpImage<float> &IDepthMask,
				   vpImage<unsigned char> &IAugmented,
           float alpha)
{
  if ((I1.getHeight() != I2.getHeight()) || (I1.getWidth() != I2.getWidth()) || (I1.getHeight() != IDepthMask.getHeight()) || (I1.getWidth() != IDepthMask.getWidth()))
  {
    throw (vpException(vpException::dimensionError, "Images have not the same size"));
  }
  
  if ((I1.getHeight() != IAugmented.getHeight()) || (I1.getWidth() != IAugmented.getWidth()))
    IAugmented.resize(I1.getHeight(), I1.getWidth(), false);

  if(alpha < 0.f)
    alpha = 0.f;
  else
    if(alpha > 1.f)
      alpha = 1.f;

  int n = I1.getHeight() * I1.getWidth() ;
  int diff ;
  unsigned char *pt_I1 = I1.bitmap, *pt_I2 = I2.bitmap, *pt_IAugmented = IAugmented.bitmap;
  float *pt_IDepthMask = IDepthMask.bitmap;
  for (int b = 0; b < n ; b++, pt_I1++, pt_I2++, pt_IAugmented++, pt_IDepthMask++)
    {
      if(*pt_IDepthMask != -1.f)
      {
        *pt_IAugmented = (unsigned char)(*pt_I1*alpha + *pt_I2*(1.0f-alpha));
      }
      else
        *pt_IAugmented = *pt_I2;
    }
}

/*!
  Compute the Zero-mean Normalized Correlation between I1 and I2
  
  \param I1 : The first image.
  \param I2 : The second image.
  \return : ZNCC of I1 and I2
*/
double CImageTools::imageZNCC(vpImage<unsigned char> &I1, 
			                        vpImage<unsigned char> &I2)
{
  if ((I1.getHeight() != I2.getHeight()) || (I1.getWidth() != I2.getWidth()))
  {
    throw (vpException(vpException::dimensionError, "The two images have not the same size"));
  }

  int n = I1.getHeight() * I1.getWidth() ;
  double meanI1 = 0., meanI2 = 0.;
	unsigned char *pt_I1 = I1.bitmap, *pt_I2 = I2.bitmap;
  vpColVector vI1(n), vI2(n);
  double *pt_vI1 = vI1.data, *pt_vI2 = vI2.data;
  for (int b = 0; b < n ; b++, pt_I1++, pt_I2++, pt_vI1++, pt_vI2++)
  {
    meanI1 += *pt_I1;
    *pt_vI1 = *pt_I1;
    meanI2 += *pt_I2;
    *pt_vI2 = *pt_I2;
	}

	meanI1 /= (double)n;
  meanI2 /= (double)n;
	vpColVector vMeanI1(n), vMeanI2(n);
	vMeanI1 = meanI1;
	vMeanI2 = meanI2;
	//vI1 -= meanI1;
	vI1 -= vMeanI1;
	vI1.normalize();
	//vI2 -= meanI2;
	vI2 -= vMeanI2;
	vI2.normalize();

	return vpColVector::dotProd(vI1,vI2);
}

/*!
  Compute the Zero-mean Normalized Correlation between I1 and I2
  
  \param I1 : The first image.
  \param I2 : The second image (already center and normalized)
  \return : ZNCC of I1 and I2
*/
double CImageTools::imageZNCC(vpImage<unsigned char> &I1, 
			                        vpImage<double> &I2)
{
  if ((I1.getHeight() != I2.getHeight()) || (I1.getWidth() != I2.getWidth()))
  {
    throw (vpException(vpException::dimensionError, "The two images have not the same size"));
  }

  int n = I1.getHeight() * I1.getWidth() ;
  double meanI1 = 0., meanI2 = 0.;
	unsigned char *pt_I1 = I1.bitmap;
  vpColVector vI1(n), vI2(n);
  double *pt_vI1 = vI1.data;
  for (int b = 0; b < n ; b++, pt_I1++, pt_vI1++)
  {
    meanI1 += *pt_I1;
    *pt_vI1 = *pt_I1;
	}

	meanI1 /= (double)n;
	vpColVector vMeanI1(n);
	vMeanI1 = meanI1;
	//vI1 -= meanI1;
	vI1 -= vMeanI1;
	vI1.normalize();

  memcpy(vI2.data, I2.bitmap, n*sizeof(double));

	return vpColVector::dotProd(vI1,vI2);
}

/*!
  Compute the Zero-mean Normalized Correlation between I1 and I2
  
  \param I1 : The first image.
  \param vI : A vector of the second image (row major, already center and normalized)
  \return : ZNCC of I1 and I2
*/
double CImageTools::imageZNCC(vpImage<unsigned char> &I1, 
			                        vpColVector &vI2zn)
{
  int n = I1.getHeight() * I1.getWidth() ;
  if (vI2zn.getRows() != n)
  {
    throw (vpException(vpException::dimensionError, "The two images have not the same size"));
  }

  vpColVector vI1zn(n);
  CImageTools::imageZN(I1, vI1zn);

	return vpColVector::dotProd(vI1zn,vI2zn);
}

/*!
  Compute the ZNCC between I1 and I2
  Jan. 2020: bug corrige I1.getWidth() a la place de I2.getWidth()
  
  \param I1 : Source image 
  \param I2 : output zero mean and normalized image
*/
void CImageTools::imageZN(vpImage<unsigned char> &I1, 
			                    vpImage<double> &I2)
{
  int n = I1.getHeight() * I1.getWidth() ;
  vpColVector vI2(n);

  CImageTools::imageZN(I1, vI2);

  I2.resize(I1.getHeight(), I1.getWidth(), false);  
  memcpy(I2.bitmap, vI2.data, n*sizeof(double));
}

/*!
  Compute the ZNCC between I1 and I2
  
  \param I1 : Source image 
  \param I2 : output zero mean and normalized image
*/
void CImageTools::imageZN(vpImage<unsigned char> &I1, 
			                    vpColVector &V)
{
  int n = I1.getHeight() * I1.getWidth() ;
  double meanI1 = 0.;
	unsigned char *pt_I1 = I1.bitmap;

  V.resize(n, false);

  double *pt_V = V.data;
  for (int b = 0; b < n ; b++, pt_I1++, pt_V++)
  {
    meanI1 += *pt_I1;
    *pt_V = *pt_I1;
	}

	meanI1 /= (double)n;

	vpColVector vMeanI1(n);
	vMeanI1 = meanI1;
	//V -= meanI1;
	V -= vMeanI1;
	V.normalize();
}

/*!
  Compute the ZN of I1
  
  \param I1 : Source image 
  \param I2 : output zero mean and normalized image
  \param Idepth_Mask : mask as depth map
*/
void CImageTools::imageZN(vpImage<unsigned char> &I1, 
			                    vpImage<double> &I2,
                          vpImage<float> &Idepth_Mask,
                          bool ZN_masked_only)
{
  int n = I1.getHeight() * I1.getWidth() ;
  vpColVector vI2(n);

  CImageTools::imageZN(I1, vI2, Idepth_Mask, ZN_masked_only);

  I2.resize(I1.getHeight(), I1.getWidth(), false);  
  memcpy(I2.bitmap, vI2.data, n*sizeof(double));
}

/*!
  Compute the ZN of I1
  
  \param I1 : Source image 
  \param V : output zero mean and normalized image intensities vector
  \param Idepth_Mask : mask as depth map
*/
void CImageTools::imageZN(vpImage<unsigned char> &I1, 
			                    vpColVector &V,
                          vpImage<float> &Idepth_Mask,
                          bool ZN_masked_only)
{
  unsigned int n = I1.getHeight() * I1.getWidth(), nConsidered = 0 ;
  double meanI1 = 0., vNorm = 0.;
	unsigned char *pt_I1 = I1.bitmap;
  float *pt_Idepth_Mask = Idepth_Mask.bitmap;
  vpColVector V_masked(n);

  V.resize(n, false);

  double *pt_V = V.data; 
  double *pt_V_masked = V_masked.data;
  for (unsigned int b = 0; b < n ; b++, pt_I1++, pt_V++, pt_Idepth_Mask++)
  {
    *pt_V = *pt_I1;
    if(*pt_Idepth_Mask > 0.f)
    {
      meanI1 += *pt_I1;
      nConsidered++;
      *pt_V_masked = *pt_I1;
      pt_V_masked++;
    }
	}

	meanI1 /= (double)nConsidered;

  pt_V_masked = V_masked.data;
  double invSqrtSumSquare = 0.;
  for (unsigned int b = 0; b < nConsidered ; b++, pt_V_masked++)
  {
    *pt_V_masked -= meanI1;
    invSqrtSumSquare += (*pt_V_masked)*(*pt_V_masked);
  }
  invSqrtSumSquare = 1./sqrt(invSqrtSumSquare);

  pt_I1 = I1.bitmap;
  pt_V = V.data; 
  if(ZN_masked_only)
  {
    pt_Idepth_Mask = Idepth_Mask.bitmap;
    double *pt_V_masked = V_masked.data;
    for (unsigned int b = 0; b < n ; b++, pt_I1++, pt_V++, pt_Idepth_Mask++)
    {
      if(*pt_Idepth_Mask > 0.f)
      {
        *pt_V = (*pt_V_masked)*invSqrtSumSquare;
        pt_V_masked++;
      }
      else
        *pt_V = *pt_I1;
	  }
  }
  else
  {
    for (unsigned int b = 0; b < n ; b++, pt_I1++, pt_V++)
    {
      *pt_V = ( (*pt_I1) - meanI1 )*invSqrtSumSquare;
	  }
  }
  
}

