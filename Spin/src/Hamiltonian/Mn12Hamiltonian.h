////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of Mn_12 hamiltonian                        //
//                                                                            //
//                        last modification : 28/02/2001                      //
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


#ifndef MN12HAMILTONIAN_H
#define MN12HAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ThierryChain.h"
#include "Hamiltonian/AbstractHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class Mn12Hamiltonian : public AbstractHamiltonian
{

 protected:
  
  ThierryChain Chain;

  double J1, HalfJ1;
  double J2, HalfJ2;
  double J3, HalfJ3;
  double J3p, HalfJ3p;
  double J4, HalfJ4;

  double* SzSzContributions;

 public:

  // contructor from default datas
  //
  // chain = reference on Hilbert space of the associated system
  // j1 = coupling constant J1
  // j2 = coupling constant J2
  // j3 = coupling constant J3
  // j3p = coupling constant J3p
  // j4 = coupling constant J4
  Mn12Hamiltonian(const ThierryChain& chain, double j1, double j2, double j3, double j3p, double j4);

  // destructor
  //
  ~Mn12Hamiltonian();

  // set chain
  // 
  // chain = reference on Hilbert space of the associated system
  // return value = reference on current Hamiltonian
  Mn12Hamiltonian& SetChain(const ThierryChain& chain);

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of idinces 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
		       int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // Output Stream overload
  //
  // Str = reference on output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, Mn12Hamiltonian& H);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, Mn12Hamiltonian& H);

 private:
 
  // evaluate contribution of SziSzj terms
  //
  void EvaluateSzSzContributions();

};

#endif
