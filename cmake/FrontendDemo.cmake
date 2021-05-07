cmake_minimum_required(VERSION 3.0)

set(PLOTTER_DEMO_NAME plotter_demo)
set(CELL_DEMO_NAME cell_demo)

file(GLOB_RECURSE PLOTTER_DEMO_SOURCE "tests/plotter_demo/*.cpp")
file(GLOB_RECURSE CELL_DEMO_SOURCE "tests/cell_demo/*.cpp")

include_directories(
    ${PROJECT_SOURCE_DIR},
)

add_executable(${PLOTTER_DEMO_NAME} ${PLOTTER_DEMO_SOURCE})
target_link_libraries(${PLOTTER_DEMO_NAME} taichi_isolated_core)

add_executable(${CELL_DEMO_NAME} ${CELL_DEMO_SOURCE})
target_link_libraries(${CELL_DEMO_NAME} taichi_isolated_core)
