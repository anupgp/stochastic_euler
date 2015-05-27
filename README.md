**This ```README``` is not complete.**

Stochastic Integrator : Explicit Euler Method
=============================================

This is the README file for my C++ programs that implement a stepper to perfom stochastic integration using explicit Euler method.
The codes uses some of the packages from the BOOST library.

Included are following three header files.

boost_random.hpp
----------------
C++ This header file creates the required class for generating random numbers using the BOOST library

det_model_hh_post.hpp
---------------------
C++ header file containing the HH neuron class. This is similar to the classical HH implementaion

stoch_odeint_explicit_euler.hpp
-------------------------------
C++ header file contains the implementaion of stochastic integrator class using explicit_euler method

Example
=======
Also included is the below named ```main()``` that demonstrates how to use the program: ```test_stoch_euler.cpp```

Compilation
===========
``
g++ -std=c++11 test_stoch_euler.cpp -o test_stoch_euler_run
``
Execution
=========
``
./test_stoch_euler_run
``
License
=======
This simulator library is licensed under GNU GPLv3 which is described in LICENSE file found in home directory of this project.
