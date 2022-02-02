////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//              Copyright (C) 2001-2021 Nicolas Regnault & co                 //
//                                                                            //
//                                                                            //
//                    class of abstract gradient optimisers                   //
//                                                                            //
//                        last modification : 10/12/2021                      //
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


#ifndef ABSTRACTGRADIENTOPTIMIZER_H
#define ABSTRACTGRADIENTOPTIMIZER_H


#include "config.h"
#include "MathTools/Complex.h"
#include "Abstract1DComplexTrialFunction.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include <iostream>

class RealVector;


class AbstractOptimizerData
{
}


class AbstractGradientOptimizer
{
  
 public:

  // virtual destructor
  //
  virtual ~AbstractGradientOptimizer();

  // clone function 
  //
  // return value = clone of the function 
  virtual AbstractGradientOptimizer* Clone () = 0;

  
  // set new values of the trial coefficients (keeping the initial number of parameters)
  virtual void SetTrialParameters(RealVector& coefficients, AbstractOptimizerData* OptimizerData) = 0;

  // GetNextTrialParameters:
  // NewTrialParameters - vector returning new trial parameters
  virtual void GetNextTrialParameters(RealVector& NewTrialParameters, RealVector& currentTrialParameters, AbstractOptimizerData* OptimizerData) = 0;

  virtual void GetNextTrialParameters(ComplexVector& NewTrialParameters, ComplexVector& currentTrialParameters, AbstractOptimizerData* OptimizerData) = 0;
    
};

#endif
