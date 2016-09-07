TEMPLATE = app
TARGET = tracker

CONFIG -= debug_and_release
CONFIG += debug
#CONFIG += release
QMAKE_CXXFLAGS += -std=gnu++0x # this enables compiler support for the upcoming c++ standard in which the unordered_map is included
# QMAKE_CXXFLAGS_RELEASE += -O3

QT += core \
    gui \
    opengl
HEADERS += \
    PNGReader.hpp \
    VisualizerPNGReader.hpp
SOURCES += \
    main.cpp \
    PNGReader.cpp \
    VisualizerPNGReader.cpp
FORMS += *.ui

INCLUDEPATH += \
    ../libGeneralUtils \
    ../libGeneralUtils/cvmlib/include \
    ../libRangedataUtils \
    ../libTrackingUtils \
    ../libGui3DQt
LIBS += \
    -L../libGui3DQt/Gui3DQt \
    -L../libGeneralUtils/cvmlib/lib \
    -L../libGeneralUtils \
    -L../libRangedataUtils \
    -L../libTrackingUtils \
    -lTrackingUtils \
    -lRangedataUtils \
    -lGeneralUtils \
    -lGui3DQt \
    -lboost_program_options \
    -lboost_filesystem \
    -lboost_serialization \
    -lboost_regex \
    -lboost_system \
    -lpng \
    -lglut \
    -lGLU \
    -lann \
    -llapack \
    -lblas
