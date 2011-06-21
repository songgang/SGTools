project( SGTools )

cmake_minimum_required(VERSION 2.4)


if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#
# Find ITK
#


FIND_PACKAGE(ITK REQUIRED)

IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build ITKApps without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

link_directories(
  ${ITK_LIBRARY_PATH}
)

# Manage FFT v3 Options
#
# USE_FFTWD -- use double precision fftw
option(USE_FFTWD "Use double precision fftw if found" OFF)
# mark_as_advanced(USE_FFTWD)

#
# USE_FFTWF -- use single precision fftw
option(USE_FFTWF "Use single precision fftw if found" OFF)

# mark_as_advanced(USE_FFTWF)

# USE_SYSTEM_FFTW -- locate a readybuilt fftw installation
option(USE_SYSTEM_FFTW "Use an installed version of fftw" OFF)

# mark_as_advanced(USE_SYSTEM_FFTW)


if (USE_FFTWD OR USE_FFTWF)
	if(USE_SYSTEM_FFTW)
    	find_package( FFTW )
    	link_directories(${FFTW_LIBDIR})
	else(USE_SYSTEM_FFTW)
    	link_directories(${ITK_DIR}/fftw/lib)
    	include_directories(${ITK_DIR}/fftw/include)
	endif(USE_SYSTEM_FFTW)
endif(USE_FFTWD OR USE_FFTWF)


add_executable( ShikataLineFiltering ShikataLineFiltering.cxx )
target_link_libraries( ShikataLineFiltering ${ITK_LIBRARIES})



