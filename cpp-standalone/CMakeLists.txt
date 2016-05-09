cmake_minimum_required (VERSION 3.1.0) 
project (whatshapcpp CXX)

list(APPEND CMAKE_MODULE_PATH .) 

IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# enable_testing()

# find_package(Boost 1.58 COMPONENTS program_options iostreams unit_test_framework REQUIRED)


# SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
# SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -g -ggdb -O0 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -O3")
set(CMAKE_CXX_FLAGS_PROFILE "-Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -O3 -g -ggdb")
include_directories("../cereal-1.1.2/include")

add_subdirectory(../src ${CMAKE_CURRENT_BINARY_DIR}/bin)

# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS tests)
