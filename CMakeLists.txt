SET(PROJECT_NAME pose3D)
PROJECT(${PROJECT_NAME})

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
#set(CMAKE_VERBOSE_MAKEFILE ON)

set(BUILD_TYPE Release)
#set(BUILD_TYPE RelWithDebInfo)

set(EXECUTABLE_OUTPUT_PATH bin)
set(LIBRARY_OUTPUT_PATH lib)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


# required libraries

find_package(OpenGL REQUIRED)
find_package(FLTK REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenCV REQUIRED)


include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${OpenGL_INCLUDE_DIR})
include_directories(${FLTK_INCLUDE_DIR})  
include_directories(${EIGEN3_INCLUDE_DIR}) 



set(CMAKE_CXX_FLAGS
   "-O2 -g -std=c++0x"
#   "${SSE_FLAGS} -O3 -g -std=c++0x -fno-omit-frame-pointer"
)

add_subdirectory(src)
# Sources files

#set_property( TARGET dso APPEND_STRING PROPERTY COMPILE_FLAGS -Wall )
add_executable(pose3d_main ${PROJECT_SOURCE_DIR}/src/mocapPlayer.cpp)

target_link_libraries(pose3d_main pose3d)
target_link_libraries(pose3d_main ${FLTK_LIBRARIES})
target_link_libraries(pose3d_main ${OPENGL_LIBRARIES})
target_link_libraries(pose3d_main ${OpenCV_LIBS})

