cmake_minimum_required(VERSION 3.1.0)

project(TrackingUtils)

# As moc files are generated in the binary dir, tell CMake to look for includes there
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)
# Create code from a list of Qt designer ui files
set(CMAKE_AUTOUIC ON)

#set(AUTOGEN_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}")
#^^^^this seems not to work. files are generated within:
set(AUTOGEN_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_autogen")
#message(STATUS "CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}")
#message(STATUS "AUTOGEN_BUILD_DIR=${AUTOGEN_BUILD_DIR}")

find_package(Boost REQUIRED COMPONENTS filesystem thread)
find_package(Qt5 REQUIRED COMPONENTS Widgets)
find_library(ann_LIBRARIES NAMES ann)

add_library(${PROJECT_NAME} 
    AlgoFrameFeatures.cpp
    AlgoFrameFeatures.hpp
    AlgoRegistration.cpp
    AlgoRegistration.hpp
    AlgoSegmentation.cpp
    AlgoSegmentation.hpp
    AlgoTrackGeneration.cpp
    AlgoTrackGeneration.hpp
    AlgoVisualization.cpp
    AlgoVisualization.hpp
    ICPBesl.hpp
    ICPBesl.tcc
    ICPChen.hpp
    ICPChen.tcc
    ICP.cpp
    ICP.hpp
    ICPLinearized.hpp
    ICPLinearized.tcc
    ICPPlanePlane.hpp
    ICPPlanePlane.tcc
    ICPSwitch.hpp
    ICPSwitch.tcc
    KalmanFilter.cpp
    KalmanFilter.hpp
    LidarFrameBuffer.cpp
    LidarFrameBuffer.hpp
    LidarFrame.cpp
    LidarFrame.hpp
    ObstacleTracking.cpp
    ObstacleTracking.hpp
    ParameterHeap.cpp
    ParameterHeap.hpp
    PointCloudTrack.cpp
    PointCloudTrack.hpp
    TrackingUtils.cpp
    TrackingUtils.hpp
    VisualizerObstacleTracking.cpp
    VisualizerObstacleTracking.hpp
    VisualizerObstacleTracking.ui
)

qt5_use_modules(${PROJECT_NAME} Widgets)

target_include_directories(${PROJECT_NAME} PUBLIC
  ${CMAKE_CURRENT_LIST_DIR}
  ${AUTOGEN_BUILD_DIR}/include # generated ui_ files need to be exported
)

target_link_libraries(${PROJECT_NAME} 
    ${Boost_LIBRARIES}
    ${ann_LIBRARIES}
    Gui3DQt
    GeneralUtils
    RangedataUtils
)

target_compile_options(${PROJECT_NAME} PUBLIC
    -msse
    -msse2 
    -msse3
    #-O2 -ggdb3 # full debug information 
    #-O3 -D NDEBUG -D BOOST_DISABLE_ASSERTS #full optimization
)
