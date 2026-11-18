cmake_minimum_required(VERSION 3.16)
project(solve_ODE LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -ffp-contract=off -fno-omit-frame-pointer")

# Place runtime binaries in the source folder so the executable appears in the workspace root.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR/build}")

# Set the base Homebrew path for macOS (Apple Silicon default)
# This lets find_package locate Boost, Matplot++, and Indicators automatically
list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew")

# Find required packages
# find_package(Boost REQUIRED)
find_package(Matplot++ REQUIRED)
find_package(indicators REQUIRED) # Uncomment if you need indicators target

add_executable(solve_ODE ./src/main.cpp)

# Link Matplot++ and indicators; both expose their own include directories and dependencies.
target_link_libraries(solve_ODE PRIVATE
     Matplot++::matplot
     indicators::indicators
)
