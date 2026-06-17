# CMakeLists.txt — LAr TPC Simulation with Image Output
# ===========================================================
cmake_minimum_required(VERSION 3.16...3.27)
project(LArTPC_ImageSim VERSION 2.0 LANGUAGES CXX)

# Find Geant4
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Source files
file(GLOB_RECURSE sources "${PROJECT_SOURCE_DIR}/src/*.cc")
file(GLOB_RECURSE headers "${PROJECT_SOURCE_DIR}/include/*.hh")

# Executable
add_executable(LArTPCSim main.cc ${sources} ${headers})
target_link_libraries(LArTPCSim ${Geant4_LIBRARIES})

# Tüm macro (.mac) dosyalarını build klasörüne kopyala
file(GLOB MACRO_FILES "${PROJECT_SOURCE_DIR}/*.mac")
foreach(macro ${MACRO_FILES})
    get_filename_component(macro_name ${macro} NAME)
    configure_file(${macro} ${PROJECT_BINARY_DIR}/${macro_name} COPYONLY)
endforeach()

# Tüm Python (.py) scriptlerini build klasörüne kopyala (Daha pratik)
file(GLOB PYTHON_SCRIPTS "${PROJECT_SOURCE_DIR}/*.py")
foreach(script ${PYTHON_SCRIPTS})
    get_filename_component(script_name ${script} NAME)
    configure_file(${script} ${PROJECT_BINARY_DIR}/${script_name} COPYONLY)
endforeach()

# Release mode for batch production
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Geant4 version: ${Geant4_VERSION}")
message(STATUS "Project Path: ${PROJECT_SOURCE_DIR}")