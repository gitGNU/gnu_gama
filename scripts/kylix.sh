#!/bin/sh
#
# script for compiling Gama with Kylix3 compiler bc++

export PATH=/usr/local/kylix3/bin:$PATH
export CC=bc++
export CFLAGS="-w-8008 -w-8065 -w-8066 -w-8057 -w-8004"
export CXX=bc++
export CXXFLAGS=-w


make dep-expat-1.1
make clean
make 
