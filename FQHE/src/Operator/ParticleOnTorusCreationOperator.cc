////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class density operator                          // 
//              for particle on a torus with magentic translations            //
//                                                                            //
//                        last modification : 16/07/2015                      //
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
#include "Operator/ParticleOnTorusCreationOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

using std::cout;
using std::endl;


// constructor
//
// particle = hilbert space associated to the particles with the small number of particles
// creationIndex = index of the creation operator
ParticleOnTorusCreationOperator::ParticleOnTorusCreationOperator(ParticleOnTorus* particle, int creationIndex)
{
  this->Particle = (ParticleOnTorus*) (particle->Clone());
  this->CreationIndex = creationIndex;
}

// copy constructor
//

ParticleOnTorusCreationOperator::ParticleOnTorusCreationOperator(ParticleOnTorusCreationOperator& oper)
{
  this->Particle = (ParticleOnTorus*) (oper.Particle->Clone());  
  this->CreationIndex = oper.CreationIndex;
}

// destructor
//

ParticleOnTorusCreationOperator::~ParticleOnTorusCreationOperator()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnTorusCreationOperator::Clone ()
{
  return new ParticleOnTorusCreationOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusCreationOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnTorus*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnTorusCreationOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnTorusCreationOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix, located in the target space
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnTorusCreationOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  int Last = firstComponent + nbrComponent;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  int DimT = this->Particle->GetTargetSpace()->GetHilbertSpaceDimension();
  int NbrFluxQuanta = this->Particle->GetNbrOrbitals();
  if (V2.GetVectorDimension() != Dim || V1.GetVectorDimension() != DimT)
    {
      std::cout << "Vectors have wrong dimensions for given spaces in ParticleOnTorusCreationOperator::PartialMatrixElement"<<std::endl;
      return 0.0;
    }
  Complex Element = 0.0;
  double Coefficient = 0.0;
  int NbrTranslations;
  int Index;
  
  for (int i = (int) firstComponent; i < Last; ++i)
    {
      Coefficient=1.0;
      int targetIdx = this->Particle->Ad(i, this->CreationIndex, Coefficient);
      if (targetIdx < DimT)
	Element += Conj(V1[targetIdx]) * V2[i] * Coefficient;
    }
  return Element;
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnTorusCreationOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
											   int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  int DimT = this->Particle->GetTargetSpace()->GetHilbertSpaceDimension();
  int NbrFluxQuanta = this->Particle->GetNbrOrbitals();

  double Coefficient = 0.0;
  int Index;

  for (int i = (int) firstComponent; i < Last; ++i)
    {
      Coefficient=1.0;
      int targetIdx = this->Particle->Ad(i, this->CreationIndex, Coefficient);
      if (targetIdx < DimT)
	vDestination[targetIdx] += vSource[i] * Coefficient;
    }
  return vDestination;
}

