TEMPLATE = app
TARGET = visPointcloud

CONFIG += debug
CONFIG -= debug_and_release
QMAKE_CXXFLAGS += -std=gnu++0x # this enables compiler support for the upcoming c++ standard in which the unordered_map is included
#QMAKE_CXXFLAGS += -O3 -D NDEBUG -D BOOST_DISABLE_ASSERTS

QT += core \
    gui \
    opengl
HEADERS += ConsoleProgressBar.hpp Visualizer3DMap.hpp
SOURCES += ConsoleProgressBar.cpp  Visualizer3DMap.cpp main.cpp  
FORMS += *.ui

INCLUDEPATH += \
    ../libGeneralUtils \
    ../libGui3DQt
LIBS += \
    -L../libGui3DQt/Gui3DQt \
    -L../libGeneralUtils \
    -lGeneralUtils \
    -lGui3DQt \
    -lboost_program_options \
    -lboost_filesystem \
    -lboost_system \
    -lGLU \
    -llapack
#QMAKE_CXXFLAGS += -pg
