add_executable(mytest mytest_main.cpp common/mytest.cpp)

if (NOT MOTION_BUILD_BOOST_FROM_SOURCES)
    find_package(Boost
            COMPONENTS
            program_options
            REQUIRED)
endif ()

target_link_libraries(mytest
        MOTION::motion
        Boost::program_options
        )
