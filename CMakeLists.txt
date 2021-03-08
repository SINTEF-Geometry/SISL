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
  FILE(GLOB sisl_APPS app/*.C)
  FOREACH(app ${sisl_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} sisl)
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
  ENDFOREACH(app)
  
  FILE(GLOB sisl_APPS app/*.c)
  FOREACH(app ${sisl_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} sisl)
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY app)
  ENDFOREACH(app)
ENDIF()

IF(sisl_COMPILE_EXAMPLES OR sisl_COMPILE_VIEWER)
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
    "C:\\Users\\$ENV{USERNAME}\\Install\\include")

  find_package(GLUT REQUIRED)

  INCLUDE_DIRECTORIES(
    ${sisl_SOURCE_DIR}/viewer/include
    ${OPENGL_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}
    )

  FILE(GLOB sisl_viewer_SRCS viewer/src/*.cpp viewer/include/*.h)
  ADD_LIBRARY(sisl_viewer ${sisl_viewer_SRCS})
  IF(CMAKE_COMPILER_IS_GNUCC)
    TARGET_LINK_LIBRARIES(sisl_streaming sisl)
  ENDIF(CMAKE_COMPILER_IS_GNUCC)

  FILE(GLOB sisl_APPS3 viewer/*.cpp)
  FOREACH(app ${sisl_APPS3})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} sisl sisl_streaming sisl_viewer ${OPENGL_LIBRARIES} glut)
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY viewer)
  ENDFOREACH(app)
ENDIF()


# 'install' target

IF(WIN32)
  # Windows
  SET(CMAKE_INSTALL_PREFIX CACHE INTERNAL "")
  SET(sisl_INSTALL_PREFIX "$ENV{PROGRAMFILES}/SINTEF/sisl"
    CACHE PATH "Path to install SISL")
  INSTALL(TARGETS sisl DESTINATION ${sisl_INSTALL_PREFIX}/lib)
  INSTALL(DIRECTORY include/
    DESTINATION ${sisl_INSTALL_PREFIX}/include
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS sisl DESTINATION lib COMPONENT sisl)
  # include
  INSTALL(DIRECTORY include/
    COMPONENT sisl-dev
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)


# CPack stuff
IF(NOT GoTools_ALL_COMPILE)
  SET(CPACK_SOURCE_PACKAGE_FILE_NAME "sisl-4.6.0")
  SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;\\.swp$;\\.#;/#;.*~")
  INCLUDE(CPack)
ENDIF(NOT GoTools_ALL_COMPILE)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libsisl.so)
