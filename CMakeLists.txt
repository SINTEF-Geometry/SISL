PROJECT(sisl)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


# Include directories

INCLUDE_DIRECTORIES(
  ${sisl_SOURCE_DIR}/include
  )


# Make the sisl library

FILE(GLOB sisl_SRCS src/*.c)
ADD_LIBRARY(sisl ${sisl_SRCS})


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
