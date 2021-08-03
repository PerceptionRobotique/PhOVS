# Set SRC_subdir variable to all the files we want to parse during
# the build process. 
# Don't forget to update SRC_ALL variable if you add/remove a SRC_subdir 
# variable
#
# If you add/remove a directory, modify here

SET(SRC_PHOTOVS
  photometricVS/CImageTools.cpp
  photometricVS/CImageFilter.cpp
  photometricVS/CFeatureLuminanceOmni.cpp
  photometricVS/CFeatureLuminanceOmniFeature.cpp
  photometricVS/CFeatureLuminanceOmniFeatureDepth.cpp
  photometricVS/CFeatureLuminanceOmniFeaturesd.cpp
  photometricVS/CFeatureLuminanceOmniFeaturesdDepth.cpp
  photometricVS/CFeatureLuminanceOmniInteractionMatrix.cpp
  photometricVS/CFeatureLuminanceOmniNeighborhoods.cpp

  photometricVS/CFeatureDefocusedLuminance.cpp
)

SET(SRC_CAMERA
  camera/CCameraOmniParameters.cpp
  
  camera/CCameraThinLensParameters.cpp
)

SET (SRC_ALL
  ${SRC_PHOTOVS}
  ${SRC_CAMERA}
  )
