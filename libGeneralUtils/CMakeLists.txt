cmake_minimum_required(VERSION 3.1.0)

project(GeneralUtils)

find_package(Boost COMPONENTS serialization filesystem REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(libxml++ REQUIRED libxml++-2.6)
find_library(lapack_LIBRARIES NAMES lapack)
find_library(blas_LIBRARIES NAMES blas)

add_library(${PROJECT_NAME} 
    Classifier.cpp
    Classifier.hpp
    ClassLabel.cpp
    ClassLabelHistogram.cpp
    ClassLabelHistogram.hpp
    ClassLabel.hpp
    convert_coordinates.hpp
    CovarEllipsoidRendering.hpp
    CsvReader.cpp
    CsvReader.hpp
    deprecated.hpp
    EnumCreation.hpp
    FeatureMaxRangeNormalizer.cpp
    FeatureMaxRangeNormalizer.hpp
    FeatureNormalizer.hpp
    FeatureSet.cpp
    FeatureSet.hpp
    FeatureVector.cpp
    FeatureVector.hpp
    FeatureVector.tcc
    FMUtils.cpp
    FMUtils.hpp
    GridND.cpp
    GridND.hpp
    GridND.tcc
    HomogeneousTransformationMatrix.hpp
    kdtree2.hpp
    kogmo_time.c
    kogmo_time.h
    MatrixDefs.hpp
    NormalSampler.cpp
    NormalSampler.hpp
    ParamStructCreation.hpp
    PrincipalComponentAnalysis.cpp
    PrincipalComponentAnalysis.hpp
    serialize_unordered_map.hpp
    SvmClassifier.cpp
    SvmClassifier.hpp
    TagClassificationResult.hpp
    Tag.hpp
    TagInvalidEntries.hpp
    TagTrueClassLabel.hpp
    Timer.hpp
    TimeVal.cpp
    TimeVal.hpp
    TwoLayerGraph.hpp
    libsvm/svm.h
    libsvm/svm.cpp 
    kdtree2/src-c++/kdtree2.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_CURRENT_LIST_DIR}/cvmlib/include 
)

target_link_libraries(${PROJECT_NAME} 
    ${Boost_LIBRARIES}
    ${libxml++_LIBRARIES}
    ${lapack_LIBRARIES}
    ${blas_LIBRARIES}
)

target_compile_options(${PROJECT_NAME} PUBLIC
    -msse
    -msse2 
    -msse3
)
