add_executable(mybenchmark mybenchmark_main.cpp common/mybenchmark.cpp)

if (NOT MOTION_BUILD_BOOST_FROM_SOURCES)
    find_package(Boost
            COMPONENTS
            program_options
            REQUIRED)
endif ()

target_link_libraries(mybenchmark
        MOTION::motion
        Boost::program_options
        )
