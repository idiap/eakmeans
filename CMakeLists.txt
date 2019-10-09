cmake_minimum_required(VERSION 3.6.2)

project(eakmeans)

option(STATIC_LIBRARY "Build static library instead of shared library" OFF)
if(STATIC_LIBRARY)
	set(LIB_TYPE "STATIC")
else()
	set(LIB_TYPE "SHARED")
endif()
add_subdirectory(src)

option(BUILD_MAIN "Build test executable" ON)
if(BUILD_MAIN)
	set(EXE_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
	include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")
	link_directories("${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}" "${BLAS_FOLDER}/dll")
	add_executable(${PROJECT_NAME} ${EXE_SRC})
	target_link_libraries(${PROJECT_NAME} PUBLIC ${LIB_EAKMEANS} ${BLAS})
	add_dependencies(${PROJECT_NAME} libeakmeans)
	if(WIN32)
		set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/WHOLEARCHIVE:libeakmeans")
	endif()
endif()

option(PYTHON_KMEANS "Build Python integration of eakmeans." OFF)
if(PYTHON_KMEANS)
	execute_process(
	  COMMAND python "${PROJECT_SOURCE_DIR}/setup.py" build_ext -b lib
	  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
	  RESULT_VARIABLE py_result
	)
endif()
