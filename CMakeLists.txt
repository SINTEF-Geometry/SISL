PROJECT(sisl)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


# Include directories

INCLUDE_DIRECTORIES(
  ${sisl_SOURCE_DIR}/include
  )


# Make the sisl library

FILE(GLOB sisl_SRCS src/*.c include/*.h)
ADD_LIBRARY(sisl ${sisl_SRCS})

# This line may be needed only on AMD 64 bit... @jbt
IF(NOT WIN32)
  SET_TARGET_PROPERTIES(sisl PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF(NOT WIN32)


# Apps, examples, tests, ...?

FILE(GLOB sisl_APPS app/*.C)
FOREACH(app ${sisl_APPS})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} sisl)
ENDFOREACH(app)


# Install sisl
INSTALL(TARGETS sisl DESTINATION lib)
INSTALL(DIRECTORY include/
  DESTINATION include
  FILES_MATCHING PATTERN "*.h"
  PATTERN ".svn" EXCLUDE
)
