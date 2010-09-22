PROJECT(sisl)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


# Include directories

INCLUDE_DIRECTORIES(
  ${sisl_SOURCE_DIR}/include
  )


# Make the sisl library

FILE(GLOB sisl_SRCS src/*.c include/*.h)
ADD_LIBRARY(sisl ${sisl_SRCS})

# This line may be needed only on AMD 64 bit...
IF(NOT WIN32)
  SET_TARGET_PROPERTIES(sisl PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF(NOT WIN32)


# Set various compiler flags

IF(NOT WIN32)
  # Set warning level to '-Wall' on Linux
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
ELSE(NOT WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP8
  -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS")
ENDIF(NOT WIN32)


# Apps, examples, tests, ...?

FILE(GLOB sisl_APPS app/*.C)
FOREACH(app ${sisl_APPS})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} sisl)
ENDFOREACH(app)



# 'install' target

IF(WIN32)
  # Windows
  SET(CMAKE_INSTALL_PREFIX CACHE INTERNAL "")
  SET(sisl_INSTALL_PREFIX "$ENV{PROGRAMFILES}/sisl"
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
  INSTALL(TARGETS sisl DESTINATION lib)
  # include
  INSTALL(DIRECTORY include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)


# CPack stuff
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "sisl-4.5.0")
SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;\\.swp$;\\.#;/#;.*~")
INCLUDE(CPack)
