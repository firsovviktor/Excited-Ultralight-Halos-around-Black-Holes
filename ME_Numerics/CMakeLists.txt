cmake_minimum_required(VERSION 3.16)
project(Viktor_MC_integral LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -fno-omit-frame-pointer")

# Place runtime binaries in the source folder so the executable appears in the workspace root.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR/build}")

# Set the base Homebrew path for macOS (Apple Silicon default)
# This lets find_package locate Boost, Matplot++, and Indicators automatically
list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew")

# Find required packages
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
find_package(Boost REQUIRED)

add_executable(Viktor_MC_integral ./src/main.cpp)

# Link Matplot++ and indicators; both expose their own include directories and dependencies.
target_link_libraries(Viktor_MC_integral PRIVATE 
    Threads::Threads
    Boost::headers
)
