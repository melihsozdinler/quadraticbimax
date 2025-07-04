cmake_minimum_required(VERSION 3.10)
project(QuadraticBimax VERSION 1.0)

# Set C++17 as the standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add executable
add_executable(QuadraticBimax Bimax/QuadraticBimax.cpp)

# Set compiler flags
if(MSVC)
    target_compile_options(QuadraticBimax PRIVATE /W4)
else()
    target_compile_options(QuadraticBimax PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Install targets
install(TARGETS QuadraticBimax
    RUNTIME DESTINATION bin
) 