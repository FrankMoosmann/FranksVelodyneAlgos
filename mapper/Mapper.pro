TEMPLATE = app
TARGET = mapper

CONFIG += debug
CONFIG -= debug_and_release
QMAKE_CXXFLAGS += -std=gnu++0x # this enables compiler support for the upcoming c++ standard in which the unordered_map is included
#QMAKE_CXXFLAGS += -O3 -D NDEBUG -D BOOST_DISABLE_ASSERTS
QMAKE_CXXFLAGS += -O1 

QT += core \
    gui \
    opengl
HEADERS += *.hpp
SOURCES += \
    main.cpp \
    Frame.cpp \
    Map.cpp \
    Mapper.cpp \
    Visualizer3DMapper.cpp
FORMS += *.ui

INCLUDEPATH += \
    ../libGui3DQt \
    ../libGeneralUtils \
    ../libRangedataUtils \
    ../libTrackingUtils
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
    -lboost_regex \
    -lboost_program_options \
    -lboost_filesystem \
    -lboost_serialization \
    -lboost_system \
    -lGLU \
    -lpng \
    -lann \
    -llapack \
    -lblas
#QMAKE_CXXFLAGS += -pg
    
