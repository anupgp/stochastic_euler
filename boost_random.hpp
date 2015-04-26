/*
  boost_random.hpp : Random number class implimentaion header using BOOST functions

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


#ifndef BOOST_RANDOM_HPP_INCLUDED
#define BOOST_RANDOM_HPP_INCLUDED


#include <iostream>
#include <ctime>
#include <sys/types.h>
#include <unistd.h>

#include <boost/random/mersenne_twister.hpp> /* Header file for mt19937 engine */
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

typedef boost::random::mt19937 base_generator_type;
typedef boost::random::uniform_int_distribution<> uni_int_dist_type;
typedef boost::random::uniform_real_distribution<> uni_real_dist_type;
typedef boost::random::normal_distribution<double> norm_dist_type;

typedef boost::random::variate_generator<base_generator_type&, uni_int_dist_type> uni_int_var_gen_type;
typedef boost::random::variate_generator<base_generator_type&, uni_real_dist_type> uni_real_var_gen_type;
typedef boost::random::variate_generator<base_generator_type&, norm_dist_type> norm_var_gen_type;

class BOOST_RNG
{
private:
  long long seed;
  base_generator_type gen;
  uni_int_dist_type uni_int_dist;
  uni_real_dist_type uni_real_dist;
  norm_dist_type norm_dist;
  
  uni_int_var_gen_type rand_uni_int;
  norm_var_gen_type rand_norm;
public:
  BOOST_RNG(long long seed_ = std::time(0) +  (long long) getpid()):seed(seed_), gen(seed), 
								    uni_int_dist(0,1),rand_uni_int(gen,uni_int_dist),
								    norm_dist(0.0,1.0),rand_norm(gen,norm_dist) {};
  int uni_int(void);
  double norm_double(void);
};

// generates uniform distribution (int)
int BOOST_RNG::uni_int(void)
{
  return (rand_uni_int());
}

// generates normal distribution ( double )
double BOOST_RNG::norm_double(void)
{
  return (rand_norm());
}

#endif	// BOOST_RANDOM_HPP_INCLUDED
