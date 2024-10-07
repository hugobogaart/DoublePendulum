# taken from https://github.com/SasLuca/raylib-cmake-template

cmake_minimum_required(VERSION 3.29)
project(DoublePendulum)
# project(my_raylib_game C)
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 20)

add_compile_options()
add_link_options()


# Adding Raylib
include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE) # don't build the supplied examples
set(BUILD_GAMES    OFF CACHE BOOL "" FORCE) # don't build the supplied example games

FetchContent_Declare(
        raylib
        GIT_REPOSITORY "https://github.com/raysan5/raylib.git"
        GIT_TAG "master"
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
)

FetchContent_MakeAvailable(raylib)

add_executable(${PROJECT_NAME}
        src/main.cc
        src/smath.h
        src/Simulation.h
)

target_link_libraries(${PROJECT_NAME} PRIVATE raylib)

# Setting ASSETS_PATH
# target_compile_definitions(${PROJECT_NAME} PUBLIC ASSETS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/assets/") # Set the asset path macro to the absolute path on the dev machine
# arget_compile_definitions(${PROJECT_NAME} PUBLIC ASSETS_PATH="./assets") # Set the asset path macro in release mode to a relative path that assumes the assets folder is in the same directory as the game executable
