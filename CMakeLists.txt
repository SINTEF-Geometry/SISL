PROJECT(sisl)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Organize the project in folders (VS only?)
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)

# Include directories

INCLUDE_DIRECTORIES(
  ${sisl_SOURCE_DIR}/include
  )


# Make the sisl library

FILE(GLOB sisl_SRCS src/*.c include/*.h)
ADD_LIBRARY(sisl ${sisl_SRCS})
IF(CMAKE_COMPILER_IS_GNUCC)
    TARGET_LINK_LIBRARIES(sisl m)
ENDIF(CMAKE_COMPILER_IS_GNUCC)
SET_PROPERTY(TARGET sisl
  PROPERTY FOLDER "sisl/Libs")

# This line may be needed only on AMD 64 bit...
#IF(NOT WIN32)
#  SET_TARGET_PROPERTIES(sisl PROPERTIES COMPILE_FLAGS "-fPIC")
#ENDIF(NOT WIN32)


# Set various compiler flags
IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-but-set-variable -fPIC")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-unused-but-set-variable -fPIC")
ENDIF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
IF(WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP8 -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP8 -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS")
ENDIF(WIN32)


# Apps, examples, tests, ...?

option(sisl_COMPILE_APPS "Compile apps?" OFF)

option(sisl_COMPILE_EXAMPLES "Compile examples?" OFF)

option(sisl_COMPILE_VIEWER "Compile viewer?" OFF)

IF(sisl_COMPILE_APPS)
  FILE(GLOB sisl_APPS app/*.[cC])
  FOREACH(app ${sisl_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} sisl sisl_streaming)
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
  ENDFOREACH(app)
ENDIF()

IF(sisl_COMPILE_EXAMPLES OR sisl_COMPILE_VIEWER OR sisl_COMPILE_APPS)
  INCLUDE_DIRECTORIES(
    ${sisl_SOURCE_DIR}/streaming/include
    )

  FILE(GLOB sisl_go_SRCS streaming/src/*.cpp streaming/include/*.h)
  ADD_LIBRARY(sisl_streaming ${sisl_go_SRCS})
  # IF(CMAKE_COMPILER_IS_GNUCC)
  #     TARGET_LINK_LIBRARIES(sisl_go m)
  # ENDIF(CMAKE_COMPILER_IS_GNUCC)
ENDIF()

IF(sisl_COMPILE_EXAMPLES)
  FILE(GLOB sisl_APPS2 examples/*.cpp)
  FOREACH(app ${sisl_APPS2})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} sisl sisl_streaming)
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY examples)
  ENDFOREACH(app)
ENDIF()

IF(sisl_COMPILE_VIEWER)

  find_package(OpenGL REQUIRED)

  find_path(OPENGL_INCLUDE_DIR
            gl.h
            PATH_PREFIXES GL
            PATHS
            /usr/local/include
            /usr/include
            "$ENV{HOME}\\Install\\include" # For some reason HOME is not found in VS 2019.
            "C:\\Users\\$ENV{USERNAME}\\Install\\include"
	    "C:\\local\\include"
            )

  if(WIN32)

    find_path(GLUT_INCLUDE_DIR
              glut.h
              PATH_PREFIXES GL
              PATHS
              "$ENV{PROGRAMW6432}/Microsoft SDKs/Windows/v6.0A/Include"
              "$ENV{HOME}\\Install\\include" # For some reason HOME is not found in VS 2019.
              "C:\\Users\\$ENV{USERNAME}\\Install\\include"
	      "C:\\local\\include"
              )
    if(${MSVC_VERSION} EQUAL 1900)
      set(MSVC_NAME "msvc2015_")
      # MESSAGE("Visual Studio 2015!")
    elseif((${MSVC_VERSION} GREATER_EQUAL 1920) AND (${MSVC_VERSION} LESS 1930))
      # MESSAGE("Visual Studio 2019!")
      set(MSVC_NAME "msvc2019_")
    elseif((${MSVC_VERSION} GREATER_EQUAL 1930) AND (${MSVC_VERSION} LESS 1950))
      set(MSVC_NAME "msvc2022_")
    else()
      message("MSVC version not supported or not installed!")
    endif()
    if(CMAKE_CL_64)
      set(WIN_LIB_TYPE "64")
    else()
      set(WIN_LIB_TYPE "32")
    endif()

    find_library( GLUT_LIBRARY_RELEASE NAMES freeglut
    		  PATHS
  		  "~/Install/${MSVC_NAME}lib${WIN_LIB_TYPE}/Release"
  		  "C:/local/${MSVC_NAME}lib${WIN_LIB_TYPE}/Release"
                  ${GLUT_ROOT_PATH}/Release
                  #"~/Install/${MSVC_NAME}lib${WIN_LIB_TYPE}/Release"
                 )

    find_library( GLUT_LIBRARY_DEBUG NAMES freeglutd
    		  PATHS
                  ${GLUT_ROOT_PATH}/Debug
		  "~/Install/${MSVC_NAME}lib${WIN_LIB_TYPE}/Debug"
		  "C:/local/${MSVC_NAME}lib${WIN_LIB_TYPE}/Debug"
                  #"~/Install/${MSVC_NAME}lib${WIN_LIB_TYPE}/Debug"
                 )
    # message("GLUT_LIBRARY_RELEASE: ${GLUT_LIBRARY_RELEASE}")
    # message("GLUT_LIBRARY_DEBUG: ${GLUT_LIBRARY_DEBUG}")
    
    set(GLUT_glut_LIBRARY "")
    # message("GLUT_glut_LIBRARY: ${GLUT_glut_LIBRARY}")
    if(GLUT_LIBRARY_DEBUG)
      set(GLUT_glut_LIBRARY ${GLUT_glut_LIBRARY} debug ${GLUT_LIBRARY_DEBUG})
    endif()
    if(GLUT_LIBRARY_RELEASE)
      set(GLUT_glut_LIBRARY ${GLUT_glut_LIBRARY} optimized ${GLUT_LIBRARY_RELEASE})
    endif()
    # message("GLUT_glut_LIBRARY: ${GLUT_glut_LIBRARY}")
  else()
    find_package(GLUT REQUIRED)  
  endif()

  INCLUDE_DIRECTORIES(
    ${sisl_SOURCE_DIR}/viewer/include
    ${OPENGL_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}
    )

  # message("GLUT_INCLUDE_DIR: ${GLUT_INCLUDE_DIR}")
  # message("GLUT_INCLUDE_PATH: ${GLUT_INCLUDE_PATH}")
  # message("GLUT_glut_LIBRARY_RELEASE: ${GLUT_glut_LIBRARY_RELEASE}")
  # message("GLUT_glut_LIBRARY_DEBUG: ${GLUT_glut_LIBRARY_DEBUG}")
  # message("GLUT_glut_LIBRARY: ${GLUT_glut_LIBRARY}")
  # message("OPENGL_LIBRARIES: ${OPENGL_LIBRARIES}")
  # message("MSVC_NAME lib WIN_LIB_TYPE: ${MSVC_NAME}lib${WIN_LIB_TYPE}")

  FILE(GLOB sisl_viewer_SRCS viewer/src/*.cpp viewer/include/*.h)
  ADD_LIBRARY(sisl_viewer ${sisl_viewer_SRCS})
  IF(CMAKE_COMPILER_IS_GNUCC)
    TARGET_LINK_LIBRARIES(sisl_streaming sisl)
  ENDIF(CMAKE_COMPILER_IS_GNUCC)

  FILE(GLOB sisl_APPS3 viewer/*.cpp)
  FOREACH(app ${sisl_APPS3})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} sisl sisl_streaming sisl_viewer ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY viewer)
  ENDFOREACH(app)
ENDIF()


# 'install' target

# Set the installation prefix based on platform
IF(WIN32)
  SET(sisl_INSTALL_PREFIX "$ENV{PROGRAMFILES}/SINTEF/sisl" CACHE PATH "Path to install SISL")
ELSE()
  SET(sisl_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Path to install SISL")
ENDIF()

# === Install Targets and Export CMake Package ===
include(CMakePackageConfigHelpers)

# Install the sisl target with export
install(TARGETS sisl
    EXPORT SISLTargets
    ARCHIVE DESTINATION ${sisl_INSTALL_PREFIX}/lib
    LIBRARY DESTINATION ${sisl_INSTALL_PREFIX}/lib
    RUNTIME DESTINATION ${sisl_INSTALL_PREFIX}/bin
    INCLUDES DESTINATION ${sisl_INSTALL_PREFIX}/include
)

# Install header files
install(DIRECTORY include/
    DESTINATION ${sisl_INSTALL_PREFIX}/include
    FILES_MATCHING PATTERN "*.h"
)

# Export target info to a file
install(EXPORT SISLTargets
    FILE SISLTargets.cmake
    #NAMESPACE SISL::             # Allow target_link_libraries(MyApp SISL::sisl)
    DESTINATION ${sisl_INSTALL_PREFIX}/lib/cmake/SISL
)

# Generate and install version/config files
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/SISLConfigVersion.cmake"
    VERSION 4.6.0
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/SISLConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/SISLConfig.cmake"
    INSTALL_DESTINATION ${sisl_INSTALL_PREFIX}/lib/cmake/SISL
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SISLConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SISLConfigVersion.cmake"
    DESTINATION ${sisl_INSTALL_PREFIX}/lib/cmake/SISL
)

# CPack setup for packaging if needed
IF(NOT GoTools_ALL_COMPILE)
  SET(CPACK_SOURCE_PACKAGE_FILE_NAME "sisl-4.6.0")
  SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;\\.swp$;\\.#;/#;.*~")
  INCLUDE(CPack)
ENDIF()

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libsisl.so)
