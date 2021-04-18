cmake_minimum_required(VERSION 3.0)

set(TESTS_NAME frontend_demo)
set(DEMO_LIB_NAME demo)

file(GLOB_RECURSE FRONTEND_DEMO_SOURCE "tests/frontend_demo/*.cpp")

include_directories(
    ${PROJECT_SOURCE_DIR},
)

add_executable(${TESTS_NAME} ${FRONTEND_DEMO_SOURCE})
target_link_libraries(${TESTS_NAME} taichi_isolated_core)
add_library(${DEMO_LIB_NAME} SHARED "tests/frontend_demo/main.cpp")
target_link_libraries(${DEMO_LIB_NAME} taichi_isolated_core)