cmake_minimum_required(VERSION 3.1.0)

project(kitti_converters)

find_package(Boost REQUIRED COMPONENTS filesystem program_options)

add_executable(kitti_convert_velodyne
        convertVelodyne.cpp
        )

target_link_libraries(kitti_convert_velodyne
        ${Boost_LIBRARIES}
        RangedataUtils
        )

SET(FILES_TO_LINK_TO_BUILD_DIR
        convertKittiRawToMapping.sh
        convertKittiTrackingToMapping.sh
        convertKittiOdometry.sh
        imuheader4mapping.txt
        )

foreach(FILE ${FILES_TO_LINK_TO_BUILD_DIR})
    set(file_source "${CMAKE_CURRENT_SOURCE_DIR}/${FILE}")
    set(file_target "${CMAKE_CURRENT_BINARY_DIR}/${FILE}")
    add_custom_command(TARGET kitti_convert_velodyne POST_BUILD VERBATIM
            #COMMAND ${CMAKE_COMMAND} -E copy_if_different ${script_source} ${script_target}
            COMMAND ${CMAKE_COMMAND} -E create_symlink ${file_source} ${file_target}
            )
endforeach(FILE)
