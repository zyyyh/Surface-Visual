cmake_minimum_required(VERSION 2.8.4)
project(surface)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

SET(CMAKE_BUILD_TYPE DEBUG) #debug

if (NOT WIN32)
set(CMAKE_MODULE_PATH /usr/local/lib/cmake /usr/local/lib/x86_64-linux-gnu/cmake ${CMAKE_SOURCE_DIR}/cmake/modules)
set(CMAKE_PREFIX_PATH /usr/local/lib/cmake/glfw)
endif (NOT WIN32)


find_package (PkgConfig REQUIRED)
find_package (OpenGL REQUIRED)
if (WIN32)
  include_directories("C:/GLFW/include" "C:/glm" "C:/glew/include")
  link_directories("C:/glew/lib/Release/x64" "C:/GLFW/lib-mingw" )
else (WIN32)
  find_package (glfw3 REQUIRED)
  find_package (GLM REQUIRED)
  find_package (GLEW REQUIRED STATIC)
endif (WIN32)
if (APPLE)
  include_directories(/usr/local/include)
  find_library(COCOA_LIBRARY Cocoa REQUIRED)
  find_library(IOKIT_LIBRARY IOKit REQUIRED)
  find_library(COREVID_LIBRARY CoreVideo REQUIRED)
endif (APPLE)

find_package (MPFR REQUIRED)

include_directories( ${CMAKE_SOURCE_DIR}/deps )
include_directories( viz3d )


set(GLAD "${CMAKE_SOURCE_DIR}/deps/glad/glad.h"
         "${CMAKE_SOURCE_DIR}/deps/glad/glad.c")

file(GLOB SOURCE_FILES
 "viz3d/*.h"
 "viz3d/*.cc" 
)

## "surface" is the name of the executable
add_executable(surface ${SOURCE_FILES} ${GLAD})
  target_link_libraries (surface
          glfw
          ${OPENGL_LIBRARIES}
          ${GLEW_LIBRARY}
          ${MPFR_LIBRARIES}
          ${COCOA_LIBRARY} ${COREVID_LIBRARY} ${IOKIT_LIBRARY})

if (WIN32)
  target_link_libraries (surface
          ${OPENGL_LIBRARIES} glfw3 glew32s glu32 opengl32)
endif (WIN32)

#if (UNIX)
#  target_link_libraries (surface
#          ${OPENGL_LIBRARIES}
#          ${GLFW3_LIBRARY}
#          ${GLEW_LIBRARY}
#          Xxf86vm pthread Xrandr Xinerama Xi Xcursor)
#endif (UNIX)
