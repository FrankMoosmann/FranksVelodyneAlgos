cmake_minimum_required(VERSION 3.1.0)

project(visTracks)

# As moc files are generated in the binary dir, tell CMake to look for includes there
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)
# Create code from a list of Qt designer ui files
set(CMAKE_AUTOUIC ON)

find_package(Boost COMPONENTS filesystem program_options regex REQUIRED)
find_package(Qt5 REQUIRED COMPONENTS Widgets Core OpenGL)

include_directories(
  ${CMAKE_CURRENT_LIST_DIR}
)

add_executable(${PROJECT_NAME} 
    FrameData.cpp
    FrameData.hpp
    main.cpp
    TrackData.cpp
    TrackData.hpp
    VisualizerTracks.cpp
    VisualizerTracks.hpp
    VisualizerTracks.ui
)

qt5_use_modules(${PROJECT_NAME} Widgets)

target_link_libraries(${PROJECT_NAME} 
    ${Boost_LIBRARIES}
    TrackingUtils
    RangedataUtils
    GeneralUtils
    Gui3DQt
)

