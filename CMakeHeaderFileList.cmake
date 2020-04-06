# Set HEADER_subdir variable to all the files we want to parse during
# the build process. 
# Don't forget to update HEADER_ALL variable if you add/remove a 
# HEADER_subdir variable
#
# If you add/remove a directory, modify here

set(HEADER_PHOTOVS
  photometricVS/CImageTools.h
  photometricVS/CImageFilter.h
  photometricVS/CFeatureLuminanceOmni.h
  )

SET(HEADER_CAMERA
  camera/CCameraOmniParameters.h
)

SET (HEADER_ALL 
  ${HEADER_PHOTOVS}
  ${HEADER_CAMERA}
  )
