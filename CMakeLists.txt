cmake_minimum_required(VERSION 3.30)
project(Dereck)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(.)

add_executable(Dereck
    DataIO.cpp
    DateUtils.h
    FDOptionPricer.cpp
    FDOptionPricer.h
    main.cpp)
