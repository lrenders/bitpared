# If the Google Test framework is found, build the unit tests
find_package(GTest)

if (GTEST_FOUND)
	include_directories(${PROJECT_SOURCE_DIR}/src)
        add_executable(unittest bitparedtest.cpp)
	target_link_libraries(unittest GTest::GTest GTest::Main)
endif (GTEST_FOUND)
