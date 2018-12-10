cmake_minimum_required(VERSION 3.1.0)

project(tracker)

# As moc files are generated in the binary dir, tell CMake to look for includes there
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)
# Create code from a list of Qt designer ui files
set(CMAKE_AUTOUIC ON)

find_package(Boost REQUIRED COMPONENTS filesystem program_options regex)
find_package(Qt5 REQUIRED COMPONENTS Widgets)


add_executable(${PROJECT_NAME}
    main.cpp
    PNGReader.cpp
    PNGReader.hpp
    VisualizerPNGReader.cpp
    VisualizerPNGReader.hpp
    VisualizerPNGReader.ui
)

qt5_use_modules(${PROJECT_NAME} Widgets)

target_link_libraries(${PROJECT_NAME}
    ${Boost_LIBRARIES}
    TrackingUtils
    RangedataUtils
    GeneralUtils
    Gui3DQt 
)

