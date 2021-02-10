cmake_minimum_required(VERSION 3.5)

project(gmp-montgomery LANGUAGES C)

add_executable(gmp-montgomery curve.c curve.h main.c options.h point.c point.h)

target_link_libraries(${PROJECT_NAME} gmp)
