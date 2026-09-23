#!/bin/sh

cmake --build build

if [ $? -eq 0 ]; then
    build/instance_generator
fi