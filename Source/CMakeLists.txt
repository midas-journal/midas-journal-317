cmake_minimum_required(VERSION 2.4)

PROJECT( DMFFDLabeledPointSetRegistration )

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(Dart)

# Set up ITK
IF(USE_ITK)
  FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
      "Cannot build without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)
ENDIF(USE_ITK)

INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake) 
IF (USE_ITK_FILE) 
  INCLUDE(${USE_ITK_FILE}) 
ENDIF(USE_ITK_FILE) 


SET(IJ_INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}/JHCT
  ${CMAKE_SOURCE_DIR}/Utilities
)

INCLUDE_DIRECTORIES(${IJ_INCLUDE_DIRS})

LINK_DIRECTORIES(${ITK_LIBRARY_PATH})

IF(BUILD_TESTING)
# non-templated class -- this should be stored in a library and linked in...
SET(SOURCES "${CMAKE_SOURCE_DIR}/itkCommandLineParser" "${CMAKE_SOURCE_DIR}/itkCommandLineOption" )

ADD_EXECUTABLE(itkDMFFDLabeledPointSetRegistrationFilterTest itkDMFFDLabeledPointSetRegistrationFilterTest.cxx ${SOURCES})
TARGET_LINK_LIBRARIES(itkDMFFDLabeledPointSetRegistrationFilterTest ITKIO ITKNumerics ITKStatistics)

ADD_EXECUTABLE(CreateWarpedGridImage CreateWarpedGridImage.cxx)
TARGET_LINK_LIBRARIES(CreateWarpedGridImage ITKIO)

ADD_EXECUTABLE(ImageCompare ImageCompare.cxx)
TARGET_LINK_LIBRARIES(ImageCompare ITKIO)
ENDIF(BUILD_TESTING)

###
#  Perform testing
###

ADD_TEST(PSR_1 ${CMAKE_BINARY_DIR}/itkDMFFDLabeledPointSetRegistrationFilterTest 2 --point-sets [${CMAKE_SOURCE_DIR}/Data/smile.nii.gz,${CMAKE_SOURCE_DIR}/Data/frown.nii.gz] --transformation [3,2x2,${CMAKE_SOURCE_DIR}/Data/smile.nii.gz] --similarity [1.1,0.95,30,0,50,1,4,10] --optimization [30x20x10,1,1] --labels [1x2x3x4x5,0.1x0.1x0.1x0.1x0.1] --output ${CMAKE_BINARY_DIR}/psr_1 --verbose)

ADD_TEST(GRID_1 CreateWarpedGridImage 2 ${CMAKE_BINARY_DIR}/psr_1Warp.nii.gz ${CMAKE_BINARY_DIR}/psr_1WarpedGrid.nii.gz )

#ADD_TEST(COMPARE_GRID_1 ImageCompare ${CMAKE_BINARY_DIR}/psr_1WarpedGrid.nii.gz ${CMAKE_SOURCE_DIR}/Data/psr_1WarpedGrid.nii.gz )


ADD_TEST(PSR_2 itkDMFFDLabeledPointSetRegistrationFilterTest 3 --point-sets [${CMAKE_SOURCE_DIR}/Data/sphere.vtk,${CMAKE_SOURCE_DIR}/Data/square.vtk] --transformation [3,4x4x4,${CMAKE_SOURCE_DIR}/Data/square.nii.gz] --labels [1,0.5,1]  --similarity [1.0,0.93,10,0,50,1,5,10] --optimization [50x50x50,1,0.8] --output ${CMAKE_BINARY_DIR}/psr_2 --verbose)
#ADD_TEST(GRID_2 CreateWarpedGridImage 3 ${CMAKE_BINARY_DIR}/psr_2Warp.nii.gz ${CMAKE_BINARY_DIR}/psr_2WarpedGrid.nii.gz)
#ADD_TEST(COMPARE_GRID_2 ImageCompare ${CMAKE_BINARY_DIR}/psr_2WarpedGrid.nii.gz ${CMAKE_SOURCE_DIR}/Data/psr_2WarpedGrid.nii.gz )

