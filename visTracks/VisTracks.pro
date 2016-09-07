TEMPLATE = app
TARGET = visTracks

CONFIG += debug
CONFIG -= debug_and_release
QMAKE_CXXFLAGS += -std=gnu++0x # this enables compiler support for the upcoming c++ standard in which the unordered_map is included
#QMAKE_CXXFLAGS += -O3 -D NDEBUG -D BOOST_DISABLE_ASSERTS

QT += core \
    gui \
    opengl
HEADERS += FrameData.hpp  TrackData.hpp  VisualizerTracks.hpp
SOURCES += FrameData.cpp  TrackData.cpp  VisualizerTracks.cpp  main.cpp
FORMS += VisualizerTracks.ui

INCLUDEPATH += \
    ../libGeneralUtils \
    ../libRangedataUtils \
    ../libTrackingUtils \
    ../libGui3DQt
LIBS += \
    -L../libGui3DQt/Gui3DQt \
    -L../libGeneralUtils \
    -L../libRangedataUtils \
    -L../libTrackingUtils \
    -lTrackingUtils \
    -lRangedataUtils \
    -lGeneralUtils \
    -lGui3DQt \
    -lboost_program_options \
    -lboost_regex \
    -lboost_filesystem \
    -lboost_system \
    -lGLU \
    -lpng \
    -lann \
    -llapack \
    -lblas
#QMAKE_CXXFLAGS += -pg
