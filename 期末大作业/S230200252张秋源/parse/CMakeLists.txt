add_library(parse argparse.cpp)
target_include_directories(parse INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})