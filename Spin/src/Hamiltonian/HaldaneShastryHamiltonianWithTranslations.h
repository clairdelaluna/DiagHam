////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of Haldane-Shasry hamiltonian with translations         //
//                                                                            //
//                        last modification : 26/05/2011                      //
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


#ifndef HALDANESHASTRYHAMILTONIANWITHTRANSLATIONS_H
#define HALDANESHASTRYHAMILTONIANWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChainWithTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class HaldaneShastryHamiltonianWithTranslations : public AbstractHamiltonian
{

 protected:
  
  AbstractSpinChainWithTranslations* Chain;

  // coupling constant between spin in the xx and yy direction 
  double J;
  // coupling constant between spin in the zz direction 
  double Jz;
  // half coupling constant between spin in the xx and yy direction 
  double HalfJ;

  // coupling coefficients between two spins, array index is the distance between the two spins
  double* JCoupling;
  
  // number of spin 
  int NbrSpin;

  // array conating all matrix diagonal elements
  double* SzSzContributions;

  //array containing all the cosinus that are needed when computing matrix elements
  double* CosinusTable;
  //array containing all the sinus that are needed when computing matrix elements
  double* SinusTable;
  //array containing all the complex phase that are needed when computing matrix elements
  Complex* ExponentialTable;

 public:

  // default constructor
  //
  HaldaneShastryHamiltonianWithTranslations();

  // constructor from default datas
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpin = number of spin
  HaldaneShastryHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin);

  // destructor
  //
  ~HaldaneShastryHamiltonianWithTranslations();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set chain
  // 
  // chain = pointer on Hilbert space of the associated system
  // return value = reference on current Hamiltonian
  HaldaneShastryHamiltonianWithTranslations& SetChain(AbstractSpinChainWithTranslations* chain);

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

  // save precalculations in a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be stored
  // return value = true if no error occurs
  bool SavePrecalculation (char* fileName);

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  virtual Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  virtual Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					     int firstComponent, int nbrComponent);
 
  // return a list of left interaction operators
  //
  // return value = list of left interaction operators
  List<Matrix*> LeftInteractionOperators();  

  // return a list of right interaction operators 
  //
  // return value = list of right interaction operators
  List<Matrix*> RightInteractionOperators();  

  // Output Stream overload
  //
  // Str = reference on output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, HaldaneShastryHamiltonianWithTranslations& H);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, HaldaneShastryHamiltonianWithTranslations& H);

 protected:
 
  // evaluate all cosinus/sinus that are needed when computing matrix elements
  //
  void EvaluateCosinusTable();

  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

};

#endif
