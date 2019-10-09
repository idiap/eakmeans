cmake_minimum_required(VERSION 3.6.2)

project(libeakmeans)

# We enable C++11
set(CMAKE_CXX_STANDARD 11)

set(TARGET_ARCHI "x64")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${TARGET_ARCHI}" CACHE PATH "Output library (lib) directory")
set(LIBEAKMEANS_OUTPUT_DIR ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

option(USEBLAS "Use external BLAS library" OFF)
if(USEBLAS)
	set(BLAS_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/../../OpenBLAS/build" CACHE PATH "Folder where to find BLAS distribution")
	set(BLAS_NAME "libopenblas" CACHE STRING "Name of BLAS implementation")
endif()

file(GLOB LIB_SRC *.cpp *.hpp *.h)
get_filename_component(main_cpp ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp ABSOLUTE)
list(REMOVE_ITEM LIB_SRC "${main_cpp}")

if(USEBLAS)
	include_directories("${BLAS_FOLDER}/include")
	find_library(BLAS ${BLAS_NAME} PATHS "${BLAS_FOLDER}" PATH_SUFFIXES lib;dll)
endif()

add_library(${PROJECT_NAME} ${LIB_TYPE} ${LIB_SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories (${PROJECT_NAME} PUBLIC {CMAKE_CURRENT_SOURCE_DIR}/src)

if(USEBLAS)
	target_link_libraries(${PROJECT_NAME} PRIVATE ${BLAS})
	target_compile_definitions(${PROJECT_NAME} PUBLIC "WITHBLAS")
endif()
