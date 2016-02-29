#!/usr/sh

emcc foo.c -o foo.js \
    --memory-init-file 0 \
    -s 'EXPORTED_FUNCTIONS=["_foo"]' \
    -s 'NO_FILESYSTEM=1' \
    -s 'ALLOW_MEMORY_GROWTH=1' \
    -s 'BUILD_AS_WORKER=1' \
    -O3 \
    -std=c99 \
    -Wall
