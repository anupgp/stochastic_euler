/* Last Revised Date: 13_DEC_2014 */
/* Lastest status: All ok    */

/*
  det_model_hh_post.hpp : implimentation of HH system (similar but not exact for the classical HH equations)  
  !!!NOTE: all units are in SI system

  Copyright (C) 2015 Anup Gopalakrishna Pillai, Suhita Nadkarni Lab, IISER Pune <anupgpillai@gmail.com>

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

#ifndef DET_MODEL_HH_POST_HPP_INCLUDED
#define DET_MODEL_HH_POST_HPP_INCLUDED

#include <vector>
#include <iostream>
#include <cmath>
#include <boost/math/constants/constants.hpp>

//UNITS: all SI units - seconds, Volts, Ampere, Meters, Simenes, Farads

//DECLARATIONS FOR HH_POST_MODEL

class DET_MODEL_HH_POST
{
private:
  const double pi = boost::math::constants::pi<double>();
  double Cm = 1.0E-02;   // F/m^2 Specific capacitance
  //----Reversal potentials for Na, K, channels in mV
  double ENa = 50.0E-03;   // Reversal potential for Sodium currents in Volts
  double EK = -85.0E-03;    // Reversal potential for Potassium currents in Volts
  double ELeak = -81E-03; // Reversal potential for Leak currents in Volts

  double GNa = 130E01;   // Maximum conductance for Sodium channels S/m^2
  double GK = 36E01;    // Maximum conductance for Potassium channels S/m^2
  double GLeak = 0.3E01;  // Maximum conductance for Leak channels S/m^2

public:
  std::vector<double> X;

  double IExt = 0.0; //extern current in Amps
  double Area = 0.0;
  double Cap = 0.0;
  double V = 0.0;      // Membrane voltage V
  /* NEURON class explicit constructor */
  explicit DET_MODEL_HH_POST( std::vector<double> X_ = {0,0,0,0},double Area_ = 8.5E-12): Area(Area_), X(X_), Cap(Cm * Area) // Area in m^2
  {
  };
  //---Equations for bouton Na Channel (Ref: Dominique Engel & Peter Jonas, 2005, Neuron)
  double alpha_Na_act(const double V);
  double beta_Na_act(const double V);
  double alpha_Na_inact(const double V);
  double beta_Na_inact(const double V);
  //---Equations for bouton K Channel (Ref: Dominique Engel & Peter Jonas, 2005, Neuron)
  double alpha_K_act(const double V);
  double beta_K_act(const double V);
  //---Declaration of the ODE solver
  template <class State >
  void operator() ( const State &x, State &dxdt , const double  t );
};
//Function declarations
//------------- Na+ channels --- implimented from classical HH (1952)
double DET_MODEL_HH_POST::alpha_Na_act(const double V)
{
  double a = 0.1E03;
  double b = 45.33E-03;
  double c = 10E-03;
  double val = 0.0;
  val = (- a * ( (V + b) / (exp( -(V + b) / c ) - 1 ) ) );
  return (val);
}
double DET_MODEL_HH_POST::beta_Na_act(const double V)
{
  double a = 4;
  double b = 70.33E-03;
  double c = 18E-03;
  double val = 0.0;
  val = ( a * ( exp( - (V + b) / c )) );
  return (val);
}
double DET_MODEL_HH_POST::alpha_Na_inact(const double V)
{
  double val = 0.0;
  val = 0.07 * exp( -( V + 70.33E-03) /20E-03 );
  return (val);
}
double DET_MODEL_HH_POST::beta_Na_inact(const double V)
{
  double val = 0.0;
  val = ( 1 / ( exp( - (V + 40.33E-03) / 10E-03) + 1) );
  return (val);
}
//--------K+ channel alpha & beta
double DET_MODEL_HH_POST::alpha_K_act(const double V)
{
  double val = 0;
  val = -0.01E03 * (V + 60.33E-03) / ( exp( -(V +  60.33E-03) / 10E-03) - 1 );
  return (val);
};
double DET_MODEL_HH_POST::beta_K_act(const double V)
{
  double val = 0.0;
  val = 0.125 * exp (- (V + 70.33E-03) / 80E-03);
  return (val);
};
//-------HH_PRE class ODE Function

template <class State>
void DET_MODEL_HH_POST::operator() ( const State &x, State &dxdt , const double t )
{

  dxdt[0] = 1E03 * (( alpha_Na_act(x[3]) * (1 - x[0]) ) - (beta_Na_act(x[3]) * x[0]) ); // Na activation
  dxdt[1] = 1E03 * ( ( alpha_Na_inact(x[3]) * (1 - x[1]) ) - (beta_Na_inact(x[3]) * x[1])) ; // Na inactivation
  dxdt[2] = 1E03 * ( ( alpha_K_act(x[3]) * (1 - x[2]) ) - (beta_K_act(x[3]) * x[2])); // K activation
  double INa =   -(Area * GNa * pow(x[0],3) * x[1] * (x[3] - ENa) );
  double IK =   -(Area * GK * pow(x[2],4) * (x[3] - EK) );
  double ILeak = -(Area * GLeak *  (x[3] - ELeak) );
  dxdt[3] =  (1/Cap) * (INa + IK + ILeak  + (IExt * Area)  );
};
#endif // DET_MODEL_HH_POST_HPP_INCLUDED
