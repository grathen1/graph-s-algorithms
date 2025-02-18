cmake_minimum_required(VERSION 3.27)
project(both)

set(CMAKE_CXX_STANDARD 17)

add_executable(both main.cpp)

set(SFML_STATIC_LIBRARIES TRUE)
set(SFML_DIR D:/Завантаження/sfml/SFML-2.6.1/lib/cmake/SFML)
find_package(SFML COMPONENTS system window graphics audio network REQUIRED)

include_directories(D:/Завантаження/sfml/SFML-2.6.1/include)
target_link_libraries(both sfml-system sfml-window sfml-graphics sfml-audio)
