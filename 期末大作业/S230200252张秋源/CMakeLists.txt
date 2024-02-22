cmake_minimum_required(VERSION 3.10)
project(tetra4)

set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 20)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

find_package(nlohmann_json 3.2.0 REQUIRED)
find_package(superlu 6.0 REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
add_subdirectory(parse)

add_executable(tetra4
    inputjson.cpp 
    resultjson.cpp 
    tetra4.cpp 
    tetra4-main.cpp
    tetra4-struct.cpp
)

target_include_directories(tetra4
    PUBLIC "${PROJECT_SOURCE_DIR}"
)


target_link_libraries(tetra4
    PUBLIC parse
    PRIVATE nlohmann_json::nlohmann_json
    PRIVATE superlu::superlu
    PRIVATE Eigen3::Eigen
)

add_executable(parsehmascii
    parsehmascii.cpp
)

target_link_libraries(parsehmascii
    PUBLIC parse
    PRIVATE nlohmann_json::nlohmann_json
)