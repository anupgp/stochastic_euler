/*
  test_stoch_euler.cpp : Sample program to demonstrate the working of stochastic Euler

  Copyright (C) 2015 Anup Gopalakrishna Pillai, Suhita Nadkarni Lab, IISER, Pune <anupgpillai@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>

#include "det_model_hh_post.hpp"
#include "stoch_odeint_explicit_euler.hpp"

/* Specifying the file within <> lets C++ search for it at -I path */
#include <boost/numeric/odeint.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/iterator.hpp>

namespace bno = boost::numeric::odeint;
namespace myode = stoch_odeint_explicit_euler;

using state_type = std::vector<double>;
using boost_rk4_stepper_type = bno::runge_kutta4<state_type>;

int main(int argc, char *argv[])
{
  // initializes the 'classical' odeint stepper from BOOST
  boost_rk4_stepper_type boost_rk4;
  state_type X_hh = {0,1,0,-80.0E-03}; // {m,h,n,V}

  // initializes an HH object using the 'DET_MODEL_HH_POST' class
  det_model_hh_post det_hhp(X_hh);

  // initializes the stochastic explicit euler stepper
  myode::stoch_odeint_explicit_euler<state_type> odeint_euler;

  // the vector s can be used to add noise to each each element of the state vector
  // '0' indictaes no noise added.
  // setting all elements of the s vector reduces the stepper to simple explicit Euler!
  std::vector<unsigned> s = {0,0,0,1}; // {m,h,n,V}

  // initialization of observer object : outputs to stdout
  null_observer my_observer;

  double t = 0.0;
  const double dt = 15E-06;
  const double tmax = 100E-03;

  // integrate_const can be used to compute all the time steps at once. Use an oberver of your choise, defaults to the null_oberver which outputs to stdout
  // integrate_const(odeint_euler,det_hhp, det_hhp.X,t,tmax,dt,s,my_observer);
  while (t <= tmax) {
    odeint_euler.do_step( det_hhp,det_hhp.X,t,dt,s);
    // boost_rk4.do_step( det_hhp,det_hhp.X, t, dt); // alternatively

    if((t > 50E-03) && (t < (54E-03))) {
      det_hhp.IExt = 20E-02;
    }
    else {
      det_hhp.IExt = 0.0;
    }

    std::cout << t;
    for(unsigned i = 0; i < det_hhp.X.size(); ++i) {
      std::cout << " " << det_hhp.X[i];
    }
    std::cout << '\n';

    t = t + dt;
  }
}
