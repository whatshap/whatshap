cmake_minimum_required (VERSION 2.6) 

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
