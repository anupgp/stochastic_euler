/*
  stoch_odeint_explicit_euler.hpp : stochastic Euler implimentaion header

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


#ifndef STOCH_ODEINT_EXPLICIT_EULER_HPP_INCLUDED
#define STOCH_ODEINT_EXPLICIT_EULER_HPP_INCLUDED

#include <iostream>
#include <cmath>
#include <vector>
#include <ctime>

#include "boost_random.hpp"

namespace stoch_odeint_explicit_euler {

  template <class state_type>  
  class STOCH_ODEINT_EXPLICIT_EULER
  {
  private:
    long long seed;
    BOOST_RNG rand;
  public:
    
    STOCH_ODEINT_EXPLICIT_EULER(long long seed_ = std::time(0) +  (long long) getpid() ):seed(seed_), rand(seed) { };
    template < class system_type >
    void do_step(system_type &system, state_type &x, const double t, const double dt,  std::vector<unsigned> &s); 
  };
  template <class state_type> 
  template <class system_type>
  void STOCH_ODEINT_EXPLICIT_EULER<state_type>::do_step(system_type &system, state_type &x, const double t, const double dt, std::vector<unsigned> &s )
  {
    state_type d (x.size(),0.0);
    system(x,d,t); 		// obtain the derivative at x
    typename state_type::iterator xi = x.begin();
    typename state_type::iterator di = d.begin();
    std::vector<unsigned>::iterator si = s.begin();
    for(int i=0; i<x.size(); i++) {
      *(xi) = (*xi) + (dt *  (*di))  +  ( pow(dt,0.5) * rand.norm_double() * (*si)  * (*xi))  ; // Euler-Murayama algorithm
      ++xi;
      ++di;
      ++si;
    }
  };
//  End of stoch_odeint_explicit_euler namespace
}

// An observer that writes the time and state vector into stdout 
struct null_observer 
{
  template < class state_type>
  void operator()(const state_type &x, double t) const
  {
    std::cout << t << " " << x << std::endl;  
  }
};

// integrate function with observer object
template <class stepper_type, class system_type, class state_type, class observer_type>
void integrate_const(stepper_type stepper, system_type system, state_type &x, const double tstart, const double tend, const double dt, std::vector<unsigned> &s, observer_type observer)
{
  for(double t = tstart; t <= tend; t+=dt){
    stepper.do_step(system,x,t,dt,s);
    observer(x,t);
  }
};

// integrate function without observer object
template <class stepper_type, class system_type, class state_type>
void integrate_const(stepper_type stepper, system_type system, state_type &x, const double tstart, const double tend, const double dt, std::vector<unsigned> &s)
{
  for(double t = tstart; t <= tend; t+=dt){
    stepper.do_step(system,x,t,dt,s);
  }
};
#endif // STOCH_ODEINT_EXPLICIT_EULER_HPP_INCLUDED
