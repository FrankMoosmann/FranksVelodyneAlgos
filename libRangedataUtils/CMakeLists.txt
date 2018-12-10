cmake_minimum_required(VERSION 3.1.0)

project(RangedataUtils)

find_package(Boost REQUIRED COMPONENTS filesystem)
find_package(Qt5 REQUIRED COMPONENTS Widgets)
find_package(PNG REQUIRED)

add_library(${PROJECT_NAME} 
    LidarImage.hpp
    LidarImageFeatures.cpp
    LidarImageFeatures.hpp
    LidarImageProjector.cpp
    LidarImageProjector.hpp
    LidarImageProjectorPNG.cpp
    LidarImageProjectorPNG.hpp
    LidarImageVisualization.cpp
    LidarImageVisualization.hpp
    LidarSegment.cpp
    LidarSegment.hpp
    PngDistanceImage.hpp
)

qt5_use_modules(${PROJECT_NAME} Widgets)

target_include_directories(${PROJECT_NAME} PUBLIC
  ${CMAKE_CURRENT_LIST_DIR}
  ${PNG_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME} 
    ${Boost_LIBRARIES}
    ${PNG_LIBRARY}
    GeneralUtils
)

target_compile_options(${PROJECT_NAME} PUBLIC
    -msse
    -msse2 
    -msse3
    #-O2 -ggdb3 # full debug information 
    #-O3 -D NDEBUG -D BOOST_DISABLE_ASSERTS #full optimization
)
