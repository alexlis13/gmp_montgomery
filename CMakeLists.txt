cmake_minimum_required(VERSION 3.5)

project(gmp-montgomery LANGUAGES C)

#Write PATH
include_directories("/home/alexis/gmp-montgomery/")

#Sources
set(SOURCE_EXE main.c)

set(SOURCE_LIB curve.c)

#Headers and libraries linking
set(HEADERS curve.h)

add_library(curve STATIC ${SOURCE_LIB} ${HEADERS})

set(SOURCE_LIB point.c)

set(HEADERS point.h)

add_library(point STATIC ${SOURCE_LIB} ${HEADERS})

add_executable(${PROJECT_NAME} "main.c")

target_link_libraries(${PROJECT_NAME} curve point gmp)

include_directories("/home/alexis/gmp-montgomery/options.h")
