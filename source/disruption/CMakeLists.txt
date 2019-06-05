add_library(disruption_api traffic_reports_api.cpp line_reports_api.cpp)
add_dependencies(disruption_api protobuf_files)

if (${CMAKE_BUILD_TESTING})
    add_subdirectory(tests)
endif()
