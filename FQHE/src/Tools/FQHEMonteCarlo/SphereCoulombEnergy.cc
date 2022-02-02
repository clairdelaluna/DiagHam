#include "SphereCoulombEnergy.h"

#include <cmath>
using std::sqrt;
using std::cout;
using std::endl;


double dsqrarg;   // shift declaration to nrutil.cc to suppress warnings if unused in module that includes nrutil.h
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)


// default constructor
SphereCoulombEnergy::SphereCoulombEnergy()
{
  this->NbrFlux=0;
}

// constructor
// nbrFlux = Number of Flux piercing sphere
// width = simple model of finite layer width for 1/sqrt(w^2+r^2)
// nbrRealObservables = optionally construct multiple observables which to measure with different weights
// nbrComplexObservables  = optionally construct multiple observables which to measure with different phases
SphereCoulombEnergy::SphereCoulombEnergy(int nbrFlux, double width, int nbrRealObservables, int nbrComplexObservables)
{
  this->Type=AbstractObservable::RealObservableT;
  this->NbrFlux = nbrFlux;
  this->WidthSqr=width*width;
  cout << "nbrFlux="<<nbrFlux<<endl;
  this->Radius = sqrt(0.5*(double)NbrFlux); // the radius is also the inverse magnetic length
  if (nbrRealObservables>0)
    {
      this->RealValues = new WeightedRealObservable*[nbrRealObservable];
      for (int i=0; i<nbrRealObservable; ++i)
	this->RealValues[i] = new WeightedRealObservable();
    }
  else
    this->RealValues=NULL;
  this->NbrRealObservables=nbrRealObservables;

  if (nbrComplexObservables>0)
    {
      this->ComplexValues = new WeightedComplexObservable*[nbrComplexObservable];
      for (int i=0; i<nbrComplexObservable; ++i)
	this->ComplexValue[i] = new WeightedComplexObservable();
    }
  else
    this->ComplexValues=NULL;
  this->NbrComplexObservables=nbrComplexObservables;
  this->NbrObservations=0;
}


// destructor
SphereCoulombEnergy::~SphereCoulombEnergy()
{
  if (NbrFlux>0)
    {
      delete Values;
    }
}

// call to make an observation
void SphereCoulombEnergy::RecordValue(double weight)
{
  ++NbrObservations;
  double E= this->GetEnergy();
  this->RealValues[0]->Observe(E,weight);
}


// call to make an observation with multiple weights
void SphereCoulombEnergy::RecordValue(double *realweights, Complex *complexweights)
{
  ++NbrObservations;
  double E= this->GetEnergy();
  for (int i=0;NbrRealObservables; ++i)
    this->RealValues[i]->Observe(E,realweights[i]);
  for (int i=0;NbrComplexObservables; ++i)
    this->ComplexValues[i]->Observe(E,complexweights[i]);
}

// calculate the energy for the current configuration
double SphereCoulombEnergy::GetEnergy()
{
  int N = this->NbrParticles;
  double E=0.0;
  if (WidthSqr>0.0)
    {
      for (int i=1;i<N;i++)
	for(int j=0;j<i;j++)
	  {
	    E+=1.0/sqrt(WidthSqr+SqrNorm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]));
	  }
    }
  else
    {
      for (int i=1;i<N;i++)
	for(int j=0;j<i;j++)
	  {
	    E+=1.0/Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	  }
    }
  return 0.5*E/Radius;
}
  

// print legend to the given stream
void SphereCoulombEnergy::PrintLegend(std::ostream &output, bool all)
{
  if (all)
    {
      output << "E\t+/-";
    }
  else
    {
      output << "E\t+/-";
    }
}

// print status to the given stream
void SphereCoulombEnergy::PrintStatus(std::ostream &output, bool all)
{
  if (NbrObservations>0)
    {
      if (all)
	{
	  output << this->Values->Average()<<"\t"<<this->Values->ErrorEstimate();
	}
      else
	{
	  int tmp=output.precision();
	  output.precision(6);
	  output << this->Values->Average()<<"\t"<<this->Values->ErrorEstimate();
	  output.precision(tmp);
	}
    }
}

// print formatted data suitable for plotting
// ouput = the target stream
void SphereCoulombEnergy::WriteDataFile(std::ostream &output)
{
  output << "#  E  \t err  ";
  for (int i=1; i<this->NbrRealObservables; ++i)
    output << "  E_"<<i<<" \t err  ";
  for (int i=0; i<this->NbrComplexObservables; ++i)
    output << "  Re(E_"<<i<<") \t  Im(E_"<<i<<") \t err  ";
  output<<endl;

  output << this->RealValues[i]->Average()
	 <<"\t"<<this->RealValues[i]->ErrorEstimate();  
 
  for (int i=1; i<this->NbrRealObservables; ++i)
    output <<"\t" << this->RealValues[i]->Average()
	 <<"\t"<<this->RealValues[i]->ErrorEstimate();  

  for (int i=0; i<this->NbrComplexObservables; ++i)
    {
      Complex TmpC=this->ComplexValues[i]->Average();
      output <<"\t" << TmpC.Re <<"\t"<<TmpC.Im
	 <<"\t"<<this->ComplexValues[i]->ErrorEstimate();  

    }
}


// access the data averaged over the run:
// for real data
void SphereCoulombEnergy::GetRealWeightedEnergies(RealVector &energies)
{
  if (energies.GetVectorLength() != this->NbrRealObservables)
    std::cerr << "Warning - requested number of energies does not match observable" << endl;

  for (int i=0; i<this->NbrRealObservables; ++i)
    energies[i]= this->RealValues[i]->Average();
}

// for real errors
void SphereCoulombEnergy::GetRealWeightedErrors(RealVector &errors)
{
  if (energies.GetVectorLength() != this->NbrRealObservables)
    std::cerr << "Warning - requested number of energies does not match observable" << endl;

  for (int i=0; i<this->NbrRealObservables; ++i)
    errors[i]= this->RealValues[i]->ErrorEstimate();
}

// for complex data
void SphereCoulombEnergy::GetComplexWeightedEnergies(ComplexVector &energies)
{
  if (energies.GetVectorLength() != this->NbrComplexObservables)
    std::cerr << "Warning - requested number of energies does not match observable" << endl;

  for (int i=0; i<this->NbrComplexObservables; ++i)
    energies[i]= this->ComplexValues[i]->Average();
}


// for real errors
void SphereCoulombEnergy::GetComplexWeightedErrors(RealVector &errors)
{
  if (energies.GetVectorLength() != this->NbrComplexObservables)
    std::cerr << "Warning - requested number of energies does not match observable" << endl;

  for (int i=0; i<this->NbrComplexObservables; ++i)
    errors[i]= this->ComplexValues[i]->ErrorEstimate();
}



// set particle collection that the observable operates on
// system = particle collection
void SphereCoulombEnergy::SetParticleCollection(AbstractParticleCollection *system)
{
  if (system->GetCollectionType()!=AbstractParticleCollection::OnSphereCollection)
    {
      cout << "Need a particle collection on the sphere for SphereCoulombEnergy"<<endl;
      exit(1);
    }
  this->System = (ParticleOnSphereCollection*) system;
  this->NbrParticles = System->GetNbrParticles();
  this->System->GetSpinorCoordinates(SpinorUCoordinates, SpinorVCoordinates);
}

// additional routines for energy observables:
// returns the total background energy
double SphereCoulombEnergy::GetTotalBackgroundEnergy()
{
  return 0.5*this->NbrParticles*this->NbrParticles/this->Radius; // not tested...
}
