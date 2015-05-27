***This README is not complete.***


STOCHASTIC INTEGRATOR : EXPLICIT EULER METHOD
=============================================

This is the README file for my C++ programs that implement a stepper to perfom stochastic integration using explicit Euler method.
The codes uses some of the packages from the BOOST library.

Included are three header files:

boost_random.hpp
----------------
C++ This header file creates the required class for generating random numbers using the BOOST library

det_model_hh_post.hpp
---------------------
C++ header file containing the HH neuron class. This is similar to the classical HH implementaion

stoch_odeint_explicit_euler.hpp
-------------------------------
C++ header file contains the implementaion of stochastic integrator class using explicit_euler method

Also included is the below named 'main' function that demonstrates how to use the program

1. test_stoch_euler.cpp 

License
=======

This simulator library is licensed under GNU GPLv3 which is described in LICENSE file found in home directory of this project.
