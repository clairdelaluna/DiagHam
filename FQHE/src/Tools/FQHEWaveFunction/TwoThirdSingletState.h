////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar M�ller                  //
//                                                                            //
//                                                                            //
//           class of Jain composite fermion wave function on sphere          //
//                      with filled (pseudo) Landau levels                    //
//                                                                            //
//                        last modification : 18/05/2007                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef THOTHIRDSINGLETSTATE_H
#define THOTHIRDSINGLETSTATE_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"
#include "Matrix/ComplexMatrix.h"
#include "Tools/FQHEWaveFunction/JainCFOnSphereOrbitals.h"

class TwoThirdSingletState: public Abstract1DComplexFunctionOnSphere
{

 protected:

  JainCFOnSphereOrbitals *OrbitalFactory;

  int NbrParticles;
  int EffectiveFlux;

  GarbageFlag Flag;
  
  // factor used to multiply each element of the Slater matrix and Cauchy Permanent
  double SlaterNorm;
  double CauchyNorm;

  ComplexMatrix Orbitals;
  ComplexMatrix *Slater;
  ComplexMatrix *CauchyPermanent;

  // For internal communication with AdaptNorm:
  Complex DeterminantValue;
  Complex PermanentValue;

  // single-particle Jastrow factors
  Complex **Jij;
  
  // if particles are very close to each other, interpolation occurs in JainCFOrbitals
  // this variable is used to pass on this value between the different subroutines
  double Interpolation;
  
 public:

  // default constructor
  //
  TwoThirdSingletState();

  // constructor
  //
  // nbrParticles = number of particles
  TwoThirdSingletState(int nbrParticles);

  // copy constructor
  //
  // function = reference on the wave function to copy
  TwoThirdSingletState(const TwoThirdSingletState& function);

  // destructor
  //
  ~TwoThirdSingletState();

  // clone function 
  //
  // return value = clone of the function 
  virtual Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  virtual Complex operator ()(RealVector& x);

  // evaluate function at a given point
  //
  // uv = ensemble of spinor variables on sphere describing point
  //      where function has to be evaluated
  //      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
  // return value = function value at (uv)
  virtual Complex CalculateFromSpinorVariables(ComplexVector& uv);


  Complex GetTestValue(RealVector& x);

  // normalize the wave-function to one for the given particle positions
  // x = point where the function has to be evaluated
  void AdaptNorm(RealVector& x);

  // normalize the wave-function over an average number of MC positions
  void AdaptAverageMCNorm(int thermalize=100, int average=250);

 private:

  // do all precalculation operations required for a new set of positions

  void EvaluateTables();

};




#endif //THOTHIRDSINGLETSTATE
