cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(mimapr CXX)

#========================================

find_package(SFML COMPONENTS system window graphics CONFIG REQUIRED)

#========================================

add_executable(
	lab
	"src/main.cpp"
)

target_link_libraries(
	lab 
	PRIVATE 
	sfml-system 
	sfml-graphics 
	sfml-window
)

target_include_directories(
	lab
	PRIVATE
	"include/"
)

add_custom_target(
	resources
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/resources ${CMAKE_CURRENT_BINARY_DIR}/resources
)

add_dependencies(
	lab 
	resources
)
