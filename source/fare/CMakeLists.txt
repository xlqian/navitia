SET(GEOREF_SRC
    fare.h
    fare.cpp
)

add_library(fare ${GEOREF_SRC})

if (${CMAKE_BUILD_TESTING})
    add_subdirectory(tests)
endif()
