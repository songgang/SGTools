project( SGTools )

cmake_minimum_required(VERSION 2.4)


if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#
# Find ITK
#
include (${CMAKE_ROOT}/Modules/FindITK.cmake)
if (USE_ITK_FILE)
  include(${USE_ITK_FILE})
endif(USE_ITK_FILE)

set( CMAKE_INCLUDE_DIRECTORIES_BEFORE ON )

# set(ITK_SOURCE_DIR /home/songgang/local/project/ITKv4/ITK)

# message( ${ITK_SOURCE_DIR} )


set(CMAKE_MODULE_PATH ${ITK_SOURCE_DIR}/CMake)
option(USE_FFTWD "Use double precision FFTW if found" OFF)
option(USE_FFTWF "Use single precision FFTW if found" OFF)
find_package( FFTW )


if(USE_FFTWD)
   link_libraries(${FFTWD_LIB})
endif(USE_FFTWD)

if(USE_FFTWF)
   link_libraries(${FFTWF_LIB})
endif(USE_FFTWF)




link_directories(
  ${ITK_LIBRARY_PATH}
#  ${FFTW_LIBDIR}
)

# set(SOURCES "${CMAKE_SOURCE_DIR}/../ExperimentalITK/Common/itkCommandLineParser" "${CMAKE_SOURCE_DIR}/../ExperimentalITK/Common/itkCommandLineOption" )

# set(UI_SOURCES "${CMAKE_SOURCE_DIR}/../ExperimentalITK/Common/antsCommandLineParser" "${CMAKE_SOURCE_DIR}/../ExperimentalITK/Common/antsCommandLineOption" )

add_executable( ShikataLineFiltering ShikataLineFiltering.cxx )
target_link_libraries( ShikataLineFiltering ${ITK_LIBRARIES})

