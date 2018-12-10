cmake_minimum_required(VERSION 3.1.0)

set (CMAKE_CXX_STANDARD 11)

add_subdirectory(libGeneralUtils)
add_subdirectory(libGui3DQt/Gui3DQt)
add_subdirectory(libRangedataUtils)
add_subdirectory(libTrackingUtils)
add_subdirectory(mapper)
add_subdirectory(tracker)
add_subdirectory(visPointCloud)
add_subdirectory(visTracks)
