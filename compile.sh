#!/bin/bash

clean=$1

if [[ $clean == "clean" ]]; then
	rm -rf bin build deps
	git submodule update --recursive --init
fi

#cmake -DCMAKE_C_COMPILER="gcc-10" -DCMAKE_CXX_COMPILER=g++-10 -DBUILD_STATIC=1 -H. -Bbuild && cmake --build build -- -j 30
cmake -DCMAKE_C_COMPILER="gcc-10" -DCMAKE_CXX_COMPILER=g++-10 -H. -Bbuild && cmake --build build -- -j 30
#cmake -H. -Bbuild && cmake --build build -- -j 5
