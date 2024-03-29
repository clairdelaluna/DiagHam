////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//      class for a basic Monte Carlo algorith for particles on a sphere      //
//                                                                            //
//                        last modification : 23/01/2008                      //
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


#include "SimpleMonteCarloAlgorithm.h"
#include "TrivialSamplingFunction.h"
#include "Options/Options.h"
#include "ParticleOnSphereCollection.h"
#include "ParticleOnDiskCollection.h"

#include <iostream>
using std::cout;
using std::endl;

// default constructor
SimpleMonteCarloAlgorithm::SimpleMonteCarloAlgorithm()
{
  this->NbrParticles=0;
  this->HavePrivateSamplingFct=false;
}

// set up for basic monte-carlo scheme
// nbrParticles = number of particles in system
// waveFunction = wavefunction to be simulated
// samplingFunction = function to be used to generate samples
// manager = pointer to option manager
// maxNbrObservables = maximum number of observables to be assigned
SimpleMonteCarloAlgorithm::SimpleMonteCarloAlgorithm(AbstractParticleCollection::Types type, int nbrParticles,
		    Abstract1DComplexFunction *waveFunction, AbstractMCSamplingFunction *samplingFunction,
								     OptionManager *manager, int maxNbrObservables, double nu)
{
  if (waveFunction==0)
    {
      cout << "Invalid wavefunction" << endl;
      exit(1);
    }
  if (samplingFunction==0)
    {
      cout << "Attention, using inefficient sampling from wavefunction itself" << endl;
      this->SamplingFunction=new TrivialSamplingFunction(waveFunction);
      HavePrivateSamplingFct=true;
    }
  else
    {
      this->SamplingFunction=samplingFunction;
      HavePrivateSamplingFct=false;
    }
  this->NbrParticles=nbrParticles;
  this->WaveFunction=waveFunction; // is returned normalized by wavefunction handler...

  this->Options=manager;  
  long Seed = Options->GetInteger("randomSeed");
  switch (type)
    {
    case AbstractParticleCollection::OnSphereCollection:
      {
	this->System=new ParticleOnSphereCollection(this->NbrParticles, Seed);
	break;
      }
    case AbstractParticleCollection::OnDiskCollection:
      {
	this->System=new ParticleOnDiskCollection(this->NbrParticles, nu, Seed);
	break;
      }
    default:
      {
	std::cerr << "Error: unknown particle collection type" << endl;
	exit(1);
      }
    }
  this->SamplingFunction->RegisterSystem(this->System);
  this->SamplingFunction->AdaptAverageMCNorm(Options->GetInteger("thermalize")); // this also relaxes the particle positions in System
  // renormalize wavefunction, as this has led to problems
  this->NormalizePsi(100);
  this->MaxNbrObservables=maxNbrObservables;
  this->Observables= new AbstractObservable*[MaxNbrObservables];
  this->Frequencies= new int[MaxNbrObservables];
  this->NbrObservables=0;
}
  
// destructor
SimpleMonteCarloAlgorithm::~SimpleMonteCarloAlgorithm()
{
  if (this->NbrParticles!=0)
    {
      delete this->System;      
      delete [] this->Frequencies;
      if (this->HavePrivateSamplingFct)
	delete SamplingFunction;
      for (int i=0; i<NbrObservables; ++i)
	delete this->Observables[i];
      delete [] this->Observables;
    }
}

// add an observable
// O = pointer to the observable to be added
// frequency = integer value indicating on which microsteps observations are being made
void SimpleMonteCarloAlgorithm::AddObservable(AbstractObservable *O, int frequency)
{
  if (this->NbrObservables<this->MaxNbrObservables)
    {
      this->Frequencies[this->NbrObservables]=frequency;
      this->Observables[this->NbrObservables++]=O;
    }
  else
    {
      AbstractObservable **newObservables = new AbstractObservable*[this->MaxNbrObservables+10];
      int *newFrequencies = new int[this->MaxNbrObservables+10];
      for (int i=0; i<NbrObservables; ++i)
	{
	  newObservables[i]=this->Observables[i];
	  newFrequencies[i]=this->Frequencies[i];
	}
      newFrequencies[this->NbrObservables]=frequency;
      newObservables[NbrObservables++]=O;
      delete [] this->Observables;
      this->Observables = newObservables;
      this->Frequencies=newFrequencies;
    }
  O->SetParticleCollection(this->System);
}

// thermalize system with a number of microsteps
// time = number of microsteps
// startFromRandom = flag indicating if we want to restart from a random configuration
void SimpleMonteCarloAlgorithm::Thermalize(int time, bool startFromRandom)
{
  if (startFromRandom)
    System->Randomize();
  this->PerformMicroSteps(time);
}

// renormalize wavefunction
// time = number of points to average
void SimpleMonteCarloAlgorithm::NormalizePsi(int time)
{
  double SumSqrPsiValues=0.0;
  for (int t=0; t<time; ++t)
    {
      PerformMicroSteps(1);
      SumSqrPsiValues+=SqrNorm((*(this->WaveFunction))(System->GetPositions()));
    }
  this->WaveFunction->Renormalize(sqrt((double)time/SumSqrPsiValues));  
}

// run simulation
// all options are included via the AddOptionGroup method
void SimpleMonteCarloAlgorithm::Simulate(ostream &Output)
{
  this->NbrAcceptedMoves = 0l;
  this->NbrAttemptedMoves = 0l;
  int NbrSteps = this->Options->GetInteger("nbr-iter");
  int DensityOfSamples = this->Options->GetInteger("sample-density");
  if (DensityOfSamples < 0)
    DensityOfSamples = this->NbrParticles;
  int NbrDisplay = this->Options->GetInteger("nbr-display");
  if (NbrDisplay==0) NbrDisplay=1;
  int DisplaySteps = NbrSteps / NbrDisplay;
  double SamplingAmplitude, Weight;
  Complex WaveFctValue, SamplingFctValue;
  Output <<"Step";
  for (int i=0; i<NbrObservables; ++i)
    {
      Output << "\t";
      Observables[i]->PrintLegend(Output);
    }
  Output << endl;
  int s=0;
  for (int d=0; d<NbrDisplay; ++d)
    {
      for (/* s */; s<(d+1)*DisplaySteps; ++s)
	{
	  this->PerformMicroSteps(DensityOfSamples);
	  SamplingFctValue = this->SamplingFunction->GetFunctionValue();
	  SamplingAmplitude = SqrNorm(SamplingFctValue);
	  //cout << "SamplingFctValue=" <<SamplingFctValue<<endl;
	  WaveFctValue = (*(this->WaveFunction))(System->GetPositions());
	  //cout << "WaveFctValue=" <<WaveFctValue<<endl;	  
	  Weight = SqrNorm(WaveFctValue)/SamplingAmplitude;
	  
	  //cout << "w="<<Weight<<" ratio="<<WaveFctValue/SamplingFctValue<<endl;
	  for (int i=0; i<NbrObservables; ++i)
	    if (s%Frequencies[i]==0) Observables[i]->RecordValue(Weight);
	}      
      Output <<(d+1)*DisplaySteps;
      for (int i=0; i<NbrObservables; ++i)
	{
	  Output << "\t";
	  Observables[i]->PrintStatus(Output);
	}
      Output << endl;
    }
  cout << "====== accepted moves: "<< (double)this->NbrAcceptedMoves / (double)this->NbrAttemptedMoves *100.0
       <<"% ======"<<endl;
}


// #define TESTING_MC

// perform a number of Monte-Carlo microsteps
// nbrSteps = number of steps
//
void SimpleMonteCarloAlgorithm::PerformMicroSteps(int nbrSteps)
{
  double acceptanceProbability;
  for (int i = 0; i < nbrSteps; ++i)
    {
#ifdef TESTING_MC
      Complex WaveFctValue = (*(this->WaveFunction))(System->GetPositions());
#endif
      System->Move();      
      ++this->NbrAttemptedMoves;
      acceptanceProbability = SamplingFunction->GetTransitionRatio();
#ifdef TESTING_MC
      Complex WaveFctValue2 = (*(this->WaveFunction))(System->GetPositions());
      if (fabs(acceptanceProbability  - SqrNorm(WaveFctValue2)/SqrNorm(WaveFctValue)) > 1e-10)
	{
	  std::cout << "Error with transition amplitudes: before = "<<SqrNorm(WaveFctValue)<<", after= "<<SqrNorm(WaveFctValue2) 
		    << ", ratio (predicted)=" <<SqrNorm(WaveFctValue2)/SqrNorm(WaveFctValue) << ", actual ratio="<<acceptanceProbability<<std::endl;
	}
      else
	std::cout << "No error: before = "<<SqrNorm(WaveFctValue)<<", after= "<<SqrNorm(WaveFctValue2) 
		  << ", ratio (predicted)=" <<SqrNorm(WaveFctValue2)/SqrNorm(WaveFctValue) << ", actual ratio="<<acceptanceProbability<<std::endl;
#endif
      if ((acceptanceProbability < 1.0) &&  (System->GetRandomNumber() > acceptanceProbability))
	{
	  System->RestoreMove();
	}
      else
	{
	  this->SamplingFunction->AcceptedMove();
	  ++this->NbrAcceptedMoves;
	}
    }
}


// add an option group containing all options related to the wave functions
//
// manager = pointer to the option manager
void SimpleMonteCarloAlgorithm::AddOptionGroup(OptionManager* manager)
{
 OptionGroup* MCGroup  = new OptionGroup ("Monte-Carlo options");
  (*(manager)) += MCGroup;  
  (*MCGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MCGroup) += new SingleIntegerOption  ('d', "sample-density", "number of Monte Carlo microsteps between iterations (-1=NbrParticles)", -1);
  (*MCGroup) += new SingleIntegerOption  ('\n',
					  "thermalize", "number of microsteps used for thermalization", 2500);  
  (*MCGroup) += new SingleIntegerOption  ('\n', "nbr-display",
					  "number of intermediate results displayed", 10);
  (*MCGroup) += new SingleIntegerOption  ('\n', "randomSeed", "random seed for internal random number generator", -1);
    
}

// write measurements of all observables to a given stream
// str = output stream
//
void SimpleMonteCarloAlgorithm::WriteObservations(ostream &str)
{
  for (int i=0; i<NbrObservables; ++i)
    {
      if (Observables[i]->IncludeInPrint())
	Observables[i]->WriteDataFile(str);
    }
}
