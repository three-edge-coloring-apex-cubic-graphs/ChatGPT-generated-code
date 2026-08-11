cmake_minimum_required(VERSION 3.20)
project(apex_verifier LANGUAGES CXX)

include(CTest)

add_library(apexlib
    src/core.cpp
    src/free_hom.cpp
    src/b2.cpp
    src/io.cpp
    src/configurations.cpp
    src/b3.cpp
    src/b4.cpp
    src/reducibility.cpp
    src/verify.cpp
)
target_compile_features(apexlib PUBLIC cxx_std_20)
target_include_directories(apexlib PUBLIC include)

function(apex_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endfunction()

apex_enable_warnings(apexlib)

add_executable(apex_verify src/main.cpp)
target_link_libraries(apex_verify PRIVATE apexlib)
apex_enable_warnings(apex_verify)

if(BUILD_TESTING)
    add_executable(apex_tests tests/test_main.cpp)
    target_link_libraries(apex_tests PRIVATE apexlib)
    target_compile_definitions(
        apex_tests PRIVATE APEX_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
    )
    apex_enable_warnings(apex_tests)
    add_test(NAME apex_tests COMMAND apex_tests)
endif()
