# Build and Install

## Requirements

This has been installed on several Linux flavors and OSX and should compile fine on Windows. It's C++ without system-dependent calls.

Install libraries in this order.

1. [MCRAND](https://github.com/afidd/mcrand) - This requires downloading a paywalled library for random number generation and relies upon the nvcc compiler from NVIDIA. The code can be modified to use Boost::Random or std::random instead.

2. Boost Libraries

3. [Semi-Markov](https://github.com/afidd/Semi-Markov) - The "testsir" branch is the one to use.

4. Gnu scientific library (GSL)

5. [HDF5 Data Model](http://www.hdfgroup.org/HDF5/release/obtain5.html)


Compiler

- Any C++11 compiler will work, so this includes g++, clang++ and msvc, among others.

Scripts

- Python is used in some. Relies on the h5py and matplotlib modules.

- Julia is used for comparison.jl, just for the individual distributions. Relatively unimportant. Relies on several modules, 

- R, used for plotting and density estimation.

## Linux and Mac Command line

Edit the Makefile to point to locations of installed libraries. There is no configure script.

## Windows or XCode

Look at the Makefile for the list of files and libraries each executable needs. Put them into a project in the IDE, and it will build. Ensure you build for -std=c++11.

