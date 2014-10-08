cmake_minimum_required (VERSION 2.6) 

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

find_package(Boost 1.36 COMPONENTS program_options REQUIRED)

add_executable (dp
  entry.h
  entry.cpp
  columnreader.h
  columnreader.cpp
  graycodes.h
  graycodes.cpp
  columnindexingiterator.h
  columnindexingiterator.cpp
  columnindexingscheme.h
  columnindexingscheme.cpp
  columncostcomputer.h
  columncostcomputer.cpp
  dptable.h
  dptable.cpp
  dp.cpp
)

target_link_libraries(dp ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
install(TARGETS dp RUNTIME DESTINATION bin)
