TEMPLATE = lib
CONFIG += staticlib
TARGET = TrackingUtils

CONFIG -= debug_and_release

CONFIG += debug
# CONFIG += release

#QMAKE_CXXFLAGS_RELEASE += -O3 -D NDEBUG -D BOOST_DISABLE_ASSERTS #full optimization
QMAKE_CXXFLAGS_RELEASE += -O2 -ggdb3 # full debug information 
QMAKE_CXXFLAGS_RELEASE += -msse -msse2 -msse3
QMAKE_CXXFLAGS += -std=gnu++0x   # this enables compiler support for the upcoming c++ standard in which the unordered_map is included

QT += core \
    gui \
    opengl

HEADERS += \
    AlgoFrameFeatures.hpp \
    AlgoRegistration.hpp \
    AlgoSegmentation.hpp \
    AlgoTrackGeneration.hpp \
    AlgoVisualization.hpp \
    ICPBesl.hpp \
    ICPBesl.tcc \
    ICPChen.hpp \
    ICPChen.tcc \
    ICP.hpp \
    ICPLinearized.hpp \
    ICPLinearized.tcc \
    ICPPlanePlane.hpp \
    ICPPlanePlane.tcc \
    ICPSwitch.hpp \
    ICPSwitch.tcc \
    KalmanFilter.hpp \
    LidarFrameBuffer.hpp \
    LidarFrame.hpp \
    ObstacleTracking.hpp \
    ParameterHeap.hpp \
    PointCloudTrack.hpp \
    TrackingUtils.hpp \
    VisualizerObstacleTracking.hpp
SOURCES += \
    AlgoFrameFeatures.cpp \
    AlgoRegistration.cpp \
    AlgoSegmentation.cpp \
    AlgoTrackGeneration.cpp \
    AlgoVisualization.cpp \
    ICP.cpp \
    KalmanFilter.cpp \
    LidarFrameBuffer.cpp \
    LidarFrame.cpp \
    ObstacleTracking.cpp \
    ParameterHeap.cpp \
    PointCloudTrack.cpp \
    TrackingUtils.cpp \
    VisualizerObstacleTracking.cpp
FORMS += \
    VisualizerObstacleTracking.ui


INCLUDEPATH += \
    ../libGui3DQt \
    ../libGeneralUtils \
    ../libRangedataUtils

LIBS +=
