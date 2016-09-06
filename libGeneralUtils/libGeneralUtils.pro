TEMPLATE = lib
CONFIG += staticlib
TARGET = GeneralUtils

CONFIG -= debug_and_release

CONFIG += debug
# CONFIG += release

#QMAKE_CXXFLAGS_RELEASE += -O3 -D NDEBUG -D BOOST_DISABLE_ASSERTS #full optimization
QMAKE_CXXFLAGS_RELEASE += -O2 -ggdb3 # full debug information 
QMAKE_CXXFLAGS_RELEASE += -msse -msse2 -msse3
QMAKE_CXXFLAGS += -std=gnu++0x   # this enables compiler support for the upcoming c++ standard in which the unordered_map is included

HEADERS += *.hpp libsvm/svm.h *.tcc
SOURCES += *.cpp kogmo_time.c libsvm/svm.cpp kdtree2/src-c++/kdtree2.cpp

INCLUDEPATH += \
    cvmlib/include 
LIBS += \
    -Lcvmlib/lib \
    -lboost_serialization \
    -lboost_filesystem \
    -lcvm \
    -lxml++-2.6
