add_executable(matmult matmult_main.cpp common/matmult.cpp)

if (NOT MOTION_BUILD_BOOST_FROM_SOURCES)
    find_package(Boost
            COMPONENTS
            program_options
            REQUIRED)
endif ()

target_link_libraries(matmult
        MOTION::motion
        Boost::program_options
        )
