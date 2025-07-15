cmake_minimum_required(VERSION 3.13)
project(BKPIR LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(SEAL REQUIRED)
find_package(OpenSSL REQUIRED)

add_executable(bkpir
    main.cpp   
    database.cpp
    mmkpirl.cpp
    bkpirl_single.cpp
    bkpirl_and.cpp
    bkpirl_or.cpp
    bkpirl_not.cpp
    utils.cpp
)

target_link_libraries(bkpir
    PRIVATE
        SEAL::seal
        OpenSSL::SSL
        OpenSSL::Crypto
)
