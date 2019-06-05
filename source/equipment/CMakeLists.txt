add_library(equipment_api equipment_api.cpp)
add_dependencies(equipment_api protobuf_files)

if (${CMAKE_BUILD_TESTING})
    add_subdirectory(tests)
endif()
