////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of Pfaffian wave function on sphere                 //
//                                                                            //
//                        last modification : 01/09/2004                      //
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


#include "config.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereWaveFunction.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"


// constructor
//
// nbrParticles = number of particles
// fermions = flag indicating whether to calculate bosonic or fermionic pfaffian

PfaffianOnSphereWaveFunction::PfaffianOnSphereWaveFunction(int nbrParticles, bool fermions)
{
  this->NbrParticles = nbrParticles;
  this->TmpPfaffian = ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->FermionFlag = fermions;
}

// copy constructor
//
// function = reference on the wave function to copy

PfaffianOnSphereWaveFunction::PfaffianOnSphereWaveFunction(const PfaffianOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->TmpPfaffian = ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->FermionFlag=function.FermionFlag;
}

// destructor
//

PfaffianOnSphereWaveFunction::~PfaffianOnSphereWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PfaffianOnSphereWaveFunction::Clone ()
{
  return new PfaffianOnSphereWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PfaffianOnSphereWaveFunction::operator ()(RealVector& x)
{
  Complex Tmp;
  Complex WaveFunction(1.0);
  double Theta;
  double Phi;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      Theta = x[i << 1];
      Phi = x[1 + (i << 1)];
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  Tmp.Re = sin(0.5 * (x[j << 1] - Theta)) * cos(0.5 * (Phi - x[1 + (j << 1)]));
	  Tmp.Im = sin(0.5 * (Theta + x[j << 1])) * sin(0.5 * (Phi - x[1 + (j << 1)]));
	  WaveFunction *= Tmp;
	  Tmp = 1.0 / Tmp;
	  this->TmpPfaffian.SetMatrixElement (i , j, Tmp);
	}
    }
  if (this->FermionFlag == true)
    WaveFunction *= WaveFunction;
  WaveFunction *= this->TmpPfaffian.Pfaffian();
  return WaveFunction;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)

Complex PfaffianOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  Complex Tmp;
  ComplexSkewSymmetricMatrix TmpPfaffian (this->NbrParticles);
  Complex WaveFunction(1.0);
  double Factor = M_PI * 0.5; 
  Complex TmpU;
  Complex TmpV; 
  for (int i = 0; i < this->NbrParticles; ++i)
    {      
      TmpU = uv[i << 1];
      TmpV = uv[(i << 1) + 1]; 
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{	  
	  Tmp = Factor * ((TmpU * uv[2*j+1]) - (TmpV * uv[2*j]));
	  WaveFunction *= Tmp;
	  Tmp = 1.0 / Tmp;
	  this->TmpPfaffian.SetMatrixElement (i , j, Tmp);
	}
    }
  if (this->FermionFlag == true)
    WaveFunction *= WaveFunction;
  WaveFunction *= this->TmpPfaffian.Pfaffian();
  return WaveFunction;
}
