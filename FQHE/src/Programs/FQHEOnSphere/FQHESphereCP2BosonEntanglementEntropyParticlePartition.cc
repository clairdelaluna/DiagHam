#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include "HilbertSpace/BosonOnCP2.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

// evaluate Hilbert space dimension for bosons
//
// nbrBosons = number of bosons
// nbrFluxQuanta = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension
// long BosonEvaluateHilbertSpaceDimension(int nbrBosons, int nbrFluxQuanta, int totalLz);

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz for bosons
//
// nbrBosons = number of bosons
// lzMax = two times momentum maximum value for a boson plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension
// long BosonShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereCP2BosonEntanglementEntropyParticlePartition" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-na", "minimum size of the particles whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-na", "maximum size of the particles whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-tza", "minimum values of Jz whose sectors has to be evaluated", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-tza", "maximum values of Jz whose sectors has to be evaluated (0 if equal to half the total system size)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-ya", "minimum values of Kz whose sectors has to be evaluated", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-ya", "maximum values of Kz whose sectors has to be evaluated (0 if equal to half the total system size)", -1);
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerate ground state");
//   (*SystemGroup) += new BooleanOption  ('\n', "compute-lvalue", "compute the L value of each reduced density matrix eigenstate");
  (*SystemGroup) += new BooleanOption  ('\n', "largest-lz", "only compute the largest block of the reduced density matrix (Lz=0 or 1/2)");
  (*SystemGroup) += new BooleanOption  ('\n', "positive-momenta", "only compute the positive jz and lz sectors");
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
//   (*SystemGroup) += new BooleanOption  ('\n', "use-svd", "use singular value decomposition instead of diagonalization to compute the entropy");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with partent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new BooleanOption ('\n', "lza-filter", "compute the eigenstates of the reduced density matrix only for a given value of lza");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "tza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Jz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "ya-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Kz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "na-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with na particles", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
//   (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHESphereCP2BosonEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if ((Manager.GetString("ground-file") != 0) && 
      (IsFile(Manager.GetString("ground-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("ground-file") << endl;
      return -1;
    }
  if ((Manager.GetString("degenerated-groundstate") != 0) && 
      (IsFile(Manager.GetString("degenerated-groundstate")) == false))
    {
      cout << "can't open file " << Manager.GetString("degenerated-groundstate") << endl;
      return -1;
    }


  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrFluxQuanta = Manager.GetInteger("nbr-flux"); 
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
//   bool ComputeLValueFlag = Manager.GetBoolean("compute-lvalue");
  bool EigenstateFlag = Manager.GetBoolean("density-eigenstate");
  bool LargestLSector = Manager.GetBoolean("largest-lz");
  bool PositiveSectors = Manager.GetBoolean("positive-momenta");
//   bool RealSpaceCut = Manager.GetBoolean("realspace-cut");
  int MinTza = Manager.GetInteger("min-tza");
  int MinYa = Manager.GetInteger("min-ya");
  int MaxTza = Manager.GetInteger("max-tza");
  int MaxYa = Manager.GetInteger("max-ya");
  int FilterTza = Manager.GetInteger("tza-eigenstate");
  int FilterYa = Manager.GetInteger("ya-eigenstate");
  int FilterNa = Manager.GetInteger("na-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
  bool ShowTimeFlag = Manager.GetBoolean("show-time");
//   bool SVDFlag = Manager.GetBoolean("use-svd");
  int* TotalTz = 0;
  int* TotalY = 0;
  bool Statistics = true;
  int NbrSpaces = 1;
  double* Weights =0;
  bool WeightFlag = false;
  BosonOnCP2** Spaces = 0;
  RealVector* GroundStates = 0;
  char** GroundStateFiles = 0;
  bool TzSymmetry = false;
  bool MinusTzSymmetry = false;
  bool TzZ3Symmetry = false;
//   int MaxLzA = Manager.GetInteger("max-lza");
//   int MinLzA = Manager.GetInteger("min-lza");
//   if ((ComputeLValueFlag == true) && (DensityMatrixFileName == 0))
//     {
//       cout << "compute-lvalue only valid when density-matrix is activated" << endl;
//       return - 1;
//     }

  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalTz = new int[1];
      TotalY = new int[1];
      Weights = new double[1];
      Weights[0] = 1.0;
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));            
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("degenerated-groundstate")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
        }
      NbrSpaces = DegeneratedFile.GetNbrLines();
      GroundStateFiles = new char* [NbrSpaces];
      TotalTz = new int[NbrSpaces];
      TotalY= new int[NbrSpaces];
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	  strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	}
      if (DegeneratedFile.GetNbrColumns() > 1)
	{
	  Weights = DegeneratedFile.GetAsDoubleArray(1);
	  WeightFlag = true;
	}
      else
	{
	  Weights = new double[NbrSpaces];
	  for (int i = 0; i < NbrSpaces; ++i)
	    Weights[i] = 1.0;
	} 
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalTz[i] = 0;
      TotalY[i] = 0;
      if (FQHEOnCP2FindSystemInfoFromVectorFileName(GroundStateFiles[i],
						       NbrParticles, NbrFluxQuanta, TotalTz[i], TotalY[i], TzSymmetry, MinusTzSymmetry, TzZ3Symmetry, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
      if (Statistics == true)
	{
	  cout << GroundStateFiles[i] << " is not a bosonic state" << endl;
	  return -1;
	}
      if ((NbrFluxQuanta + 1)*(NbrFluxQuanta + 2)/2 + NbrParticles > 64)
	{
	  cout  << "number of fermionic orbitals is too big to allow the storage of one state in one long integer" << endl;
	  return -1;
	}
      if ((TotalY[i] + 3*TotalTz[i] + 2*NbrParticles*NbrFluxQuanta) % 6 != 0)
	{
	  cout << "Y + 3Tz + 2N*Nphi should  be a multiple of 6" << endl;
	  return -1;
	}
      if ((TotalY[i] - 3*TotalTz[i] + 2*NbrParticles*NbrFluxQuanta) % 6 != 0)
	{
	  cout << "Y - 3Tz + 2N*Nphi should  be a multiple of 6" << endl;
	  return -1;
	}
      if ( TzSymmetry == true || MinusTzSymmetry == true || TzZ3Symmetry == true)
	{
	  cout << "Reduced density matrix cannot be computed on a Tz symmetrized vector - Convert to Nbody basis first" << endl;
	  return 0;
	}
      }
  

  GroundStates = new RealVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
      {
	cout << "can't open vector file " << GroundStateFiles[i] << endl;
	return -1;      
      }
  
  
  Spaces = new BosonOnCP2* [NbrSpaces];
  
  for (int i = 0; i < NbrSpaces; ++i)
    {
// #ifdef  __64_BITS__
//       if ((Norb + NbrParticles - 1) < 63)
// #else
// 	if ((Norb + NbrParticles - 1) < 31)	
// #endif
// 	  {
	Spaces[i] = new BosonOnCP2 (NbrParticles, NbrFluxQuanta, TotalTz[i], TotalY[i]);
// 	   }
// 	else
// 	  {
// 	    cout << << endl;
// 	  }
      
      if (Spaces[i]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	{
	  cout << "dimension mismatch between Hilbert space and ground state" << endl;
	  return 0;
	}
    }


  int MaxSubsystemNbrParticles = (NbrParticles >> 1) + (NbrParticles & 1);
  if (Manager.GetInteger("max-na") > 0)
    {
      MaxSubsystemNbrParticles = Manager.GetInteger("max-na");
    }
  
  int SubsystemNbrParticles = Manager.GetInteger("min-na");
  
  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "#  N    Tz    Y    lambda";
//       if (ComputeLValueFlag == true)
// 	DensityMatrixFile << " L^2 L";
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }
  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "partent");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was found in " << GroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  File.precision(14);
  cout.precision(14);
  
  double TotalTrace = 0.0;
  for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      
      int SubsystemMaxTotalR = SubsystemNbrParticles * NbrFluxQuanta;
      int SubsystemMinTotalR = 0;
      
//       if ((MinTza != -1) && (MinTza > SubsystemMinTotalR))
// 	SubsystemMinTotalTz = MinTza;
//       if ((MinJza != -1) && (MaxJza < SubsystemMaxTotalJz))
// 	SubsystemMaxTotalJz = MaxJza;
//       if ((MinKza != -1) && (MinKza > SubsystemMinTotalKz))
// 	SubsystemMinTotalKz = MinKza;
//       if ((MaxKza != -1) && (MaxKza < SubsystemMaxTotalKz))
// 	SubsystemMaxTotalJz = MaxKza;
//       if (LargestLSector == true)
// 	{
// 	  if (((LzMax * NbrParticles) & 1) == 0)
// 	    {
// 	      SubsystemTotalLz = 0;
// 	      SubsystemMaxTotalLz = 0;
// 	    }
// 	  else
// 	    {
// 	      SubsystemTotalLz = 1;
// 	      SubsystemMaxTotalLz = 1;
// 	    }
// 	}
//       if (PositiveSectors == true)
// 	{
// 	   SubsystemMinTotalY = 0;
// 	   
// 	  if (((NbrFluxQuanta * NbrParticles) & 1) == 0)
// 	    {
// 	      SubsystemMinTotalJz = 0;
// 	     }
// 	  else
// 	    {
// 	      SubsystemMinTotalJz = 1;
// 	    }
// 	}
//      
// 	
      for (int SubsystemTotalR = SubsystemMinTotalR; SubsystemTotalR <= SubsystemMaxTotalR; SubsystemTotalR += 1)
      {
	int SubsystemMaxTotalS = SubsystemNbrParticles * NbrFluxQuanta - SubsystemTotalR;
	int SubsystemMinTotalS = 0;
	for (int SubsystemTotalS = SubsystemMinTotalS; SubsystemTotalS <= SubsystemMaxTotalS ; SubsystemTotalS +=1)
	{
	  /*if (((SubsystemMaxTotalJz & 1) == ((SubsystemTotalJz + SubsystemTotalKz) & 1)) && 
	      (fabs(SubsystemTotalKz) + fabs(SubsystemTotalJz) <= SubsystemMaxTotalJz))*/	  
	  {
	    int SubsystemTotalTz = SubsystemTotalR - SubsystemTotalS;
	    int SubsystemTotalY = 3*(SubsystemTotalR + SubsystemTotalS) - 2*SubsystemNbrParticles*NbrFluxQuanta;
	    cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Tz=" << SubsystemTotalTz << " and subsystem total Y=" << SubsystemTotalY << endl;
	    timeval TotalStartingTime;
	    timeval TotalEndingTime;
	    if (ShowTimeFlag == true)
	      {
	      gettimeofday (&(TotalStartingTime), 0);
	      }
	  RealSymmetricMatrix PartialDensityMatrix;
	  RealMatrix PartialEntanglementMatrix;
	  
	  PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalTz, SubsystemTotalY, GroundStates[0], Architecture.GetArchitecture());
	  
	  if (WeightFlag == true)
	    PartialDensityMatrix *= Weights[0];
	  	 
	  for (int i = 1; i < NbrSpaces; ++i)
	    {
	      RealSymmetricMatrix TmpMatrix;
	      RealMatrix TmpEntanglementMatrix;
	      
	      TmpMatrix =  Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalTz, SubsystemTotalY,  GroundStates[i], Architecture.GetArchitecture());		      
	
	      if (WeightFlag == true)
		      TmpMatrix *= Weights[i];
	      if ( PartialDensityMatrix.GetNbrRow() != 0 ) 
		{
		PartialDensityMatrix += TmpMatrix;
		}
	      else
		{
		 PartialDensityMatrix = TmpMatrix;
		 }
	     }
	     
	  if (NbrSpaces > 1)
	    {
	       PartialDensityMatrix /= ((double) NbrSpaces);
	    }
	      
	  
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalEndingTime), 0);
	      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "reduced density matrix evaluated in " << Dt << "s" << endl;
	    }
	  if ((PartialDensityMatrix.GetNbrRow() > 1) || (PartialEntanglementMatrix.GetNbrRow() >= 1))
	    {
	      if (ShowTimeFlag == true)
		{
		  gettimeofday (&(TotalStartingTime), 0);
		}
	      RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
	      
		 
#ifdef __LAPACK__
	  if (LapackFlag == true)
	    {
	    if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles) && (FilterTza == SubsystemTotalTz ) && (FilterYa == SubsystemTotalY))
		{
		 RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(), PartialDensityMatrix.GetNbrRow(), true);
		 for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		    TmpEigenstates[i][i] = 1.0;
		 PartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
		 TmpDiag.SortMatrixDownOrder(TmpEigenstates);
		 char* TmpEigenstateName = new char[512];
		 int MaxNbrEigenstates = NbrEigenstates;
		 if (NbrEigenstates == 0)
		    MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
		 for (int i = 0; i < MaxNbrEigenstates; ++i)
		   {
		    if (TmpDiag[i] > 1e-14)
		      {
			sprintf (TmpEigenstateName,
			"bosons_cp2_density_n_%d_2s_%d_tz_%d_y_%d_na_%d_tza_%d_ya_%d.%d.vec",
			  NbrParticles, NbrFluxQuanta, TotalTz[0],TotalY[0],SubsystemNbrParticles, SubsystemTotalTz, SubsystemTotalY, i);
			  TmpEigenstates[i].WriteVector(TmpEigenstateName);
			}
		    }
		    delete[] TmpEigenstateName;
		  }
		  else
		    {
		      PartialDensityMatrix.LapackDiagonalize(TmpDiag);
		    }
		  }
// 	    else
// 		  {
// 		    if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles) && (FilterJza == SubsystemTotalJz ) && (FilterKza == SubsystemTotalKz ))
// 				{
// 				  RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
// 							    PartialDensityMatrix.GetNbrRow(), true);
// 				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
// 				    TmpEigenstates[i][i] = 1.0;
// 				  PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
// 				  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
// 				  char* TmpEigenstateName = new char[512];
// 				  int MaxNbrEigenstates = NbrEigenstates;
// 				  if (NbrEigenstates == 0)
// 				    MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
// 				  for (int i = 0; i < MaxNbrEigenstates; ++i)
// 				    {
// 				      if (TmpDiag[i] > 1e-14)
// 					{
// 					  sprintf (TmpEigenstateName,
// 			"bosons_sphere4d_density_n_%d_2s_%d_jz_%d_kz_%d_na_%d_jza_%d_kza_%d.%d.vec",
// 			  NbrParticles, NbrFluxQuanta, TotalJz[0],TotalKz[0],SubsystemNbrParticles, SubsystemTotalJz, SubsystemTotalKz, i);
// 					  TmpEigenstates[i].WriteVector(TmpEigenstateName);
// 					}
// 				    }
// 				  delete[] TmpEigenstateName;
// 				}
// 			      else
// 				{
// 				  PartialDensityMatrix.Diagonalize(TmpDiag);
// 				}
// 			    }
#else
			  PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
			
		    
		
		  if (DensityMatrixFileName != 0)
		    {
		      ofstream DensityMatrixFile;
		      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		      DensityMatrixFile.precision(14);
		      for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
			DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalTz << " " << SubsystemTotalY << " " << TmpDiag[i] << endl;
		      DensityMatrixFile.close();
		    }
		
	     
	      for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
		{
		  if (TmpDiag[i] > 1e-14)
		    {
		      EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
		      DensitySum +=TmpDiag[i];
		    }
		}
	      if (ShowTimeFlag == true)
		{
		  gettimeofday (&(TotalEndingTime), 0);
		  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		  cout << "diagonalization done in " << Dt << "s" << endl;
		}
	    }
	  else
	    if (PartialDensityMatrix.GetNbrRow() == 1)
	      {
		double TmpValue = PartialDensityMatrix(0,0);
		if (DensityMatrixFileName != 0)
		  {
		    ofstream DensityMatrixFile;
		    DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		    DensityMatrixFile.precision(14);
// 		    if (ComputeLValueFlag == false)
// 		      {
			DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalTz << " " << SubsystemTotalY << " " << TmpValue << endl;
// 		      }
// 		    else		      
// 		      {
// 			if (SubsystemNbrParticles == 1)
// 			  DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << " " << ((LzMax * (LzMax + 2)) / 4.0) << " " << (LzMax / 2.0) << endl;
// 			else
// 			  {
// 			    BosonOnSphereShort TmpDestinationHilbertSpace(SubsystemNbrParticles, SubsystemTotalLz, LzMax);
// 			    ParticleOnSphereSquareTotalMomentumOperator OperMomentum (&TmpDestinationHilbertSpace, LzMax);
// 			    RealVector TmpEigenstate(1);
// 			    TmpEigenstate[0] = 1.0;
// 			    double TmpSqrMomentum = (OperMomentum.MatrixElement(TmpEigenstate, TmpEigenstate)).Re;
// 			    double TmpMomentum = 0.0;
// 			    if (TmpSqrMomentum > 0.0)
// 			      TmpMomentum = (0.5 * (sqrt ((4.0 * TmpSqrMomentum) + 1.0) - 1.0));
// 			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << " " << TmpSqrMomentum << " " << TmpMomentum << endl;
// 			    if ((EigenstateFlag == true) && ((Manager.GetBoolean("lza-filter") == false) || (FilterLza == SubsystemTotalLz)) && ((NbrEigenstates == 0) || (NbrEigenstates > 0)))
// 			      {
// 				char* TmpEigenstateName = new char[512];
// 				RealVector TmpEigenstate(1);
// 				TmpEigenstate[0] = 1.0;
// 				sprintf (TmpEigenstateName,
// 					 "bosons_particlereduceddensity_na_%d_n_%d_2s_%d_lz_%d.0.vec",
// 					 SubsystemNbrParticles, NbrParticles, LzMax, SubsystemTotalLz);
// 				TmpEigenstate.WriteVector(TmpEigenstateName);
// 				delete[] TmpEigenstateName;
// 			      }
// 			  }
// 		      }
		    DensityMatrixFile.close();
		  }		  
		if (TmpValue > 1e-14)
		  {
		    EntanglementEntropy += TmpValue * log(TmpValue);
		    DensitySum += TmpValue;
		  }
	      }
	}
      }
       
    }
      File << SubsystemNbrParticles << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << endl;
      cout << "trace = " << DensitySum << endl;
      TotalTrace += DensitySum;
    }
  File.close();
//   if (RealSpaceCut == true)
//     cout <<"Total Trace = "<<TotalTrace<<endl;
  return 0;
}
