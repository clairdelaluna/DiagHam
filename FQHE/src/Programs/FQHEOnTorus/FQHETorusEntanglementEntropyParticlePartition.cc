#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/FermionOnTorus.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETorusEntanglementEntropyParticlePartition" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  Manager += SystemGroup;
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "kymax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-ky", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-na", "minimum size of the particles whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-na", "maximum size of the particles whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with partent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "kya-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Ky value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHETorusEntanglementEntropyParticlePartition -h" << endl;
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
  int KyMax = Manager.GetInteger("kymax"); 
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  bool EigenstateFlag = Manager.GetBoolean("density-eigenstate");
  int FilterKya = Manager.GetInteger("kya-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
  int* TotalKy = 0;
  bool Statistics = true;
  int NbrSpaces = 1;
  ParticleOnTorus** Spaces = 0;
  RealVector* GroundStates = 0;
  char** GroundStateFiles = 0;

  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalKy = new int[1];
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
       TotalKy = new int[NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	 }
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalKy[i] = 0;
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(GroundStateFiles[i],NbrParticles, KyMax, TotalKy[i], Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
     }


  GroundStates = new RealVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
      {
	cout << "can't open vector file " << GroundStateFiles[i] << endl;
	return -1;      
      }


  Spaces = new ParticleOnTorus* [NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (Statistics == true)
	{
	  Spaces[i] = new FermionOnTorus (NbrParticles, KyMax, TotalKy[i]);
	}
      else
	{
	  Spaces[i] = new BosonOnTorusShort (NbrParticles, KyMax, TotalKy[i]);
	}
      
      if (Spaces[i]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	{
	  cout << "dimension mismatch between Hilbert space and ground state" << endl;
	  return 0;
	}
    }

  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "#  N    Ky    lambda";
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
	  cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  File.precision(14);
  cout.precision(14);

  int MaxSubsystemNbrParticles = (NbrParticles >> 1) + (NbrParticles & 1);
  if (Manager.GetInteger("max-na") > 0)
    MaxSubsystemNbrParticles = Manager.GetInteger("max-na");
  int SubsystemNbrParticles = Manager.GetInteger("min-na");

  for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;

      int SubsystemMaxTotalKy = KyMax - 1;

      int SubsystemTotalKy = 0; 
      for (; SubsystemTotalKy <= SubsystemMaxTotalKy; ++SubsystemTotalKy)
	{
	  cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Ky=" << SubsystemTotalKy << endl;
	  RealSymmetricMatrix PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKy, GroundStates[0]);
	  for (int i = 1; i < NbrSpaces; ++i)
	    {
	      RealSymmetricMatrix TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalKy, GroundStates[i]);
	      PartialDensityMatrix += TmpMatrix;
	    }
	  if (NbrSpaces > 1)
	    PartialDensityMatrix /= ((double) NbrSpaces);
	  if (PartialDensityMatrix.GetNbrRow() > 1)
	    {
	      RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
	      if (LapackFlag == true)
		{
		  if ((EigenstateFlag == true) && (FilterKya == SubsystemTotalKy))
		    {
		      RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
						PartialDensityMatrix.GetNbrRow(), true);
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
			      if (Statistics == true)
				{
				  sprintf (TmpEigenstateName,
					   "fermions_torus_kysym_density_n_%d_2s_%d_ky_%d_na_%d_kya_%d.%d.vec",
					   NbrParticles, KyMax, TotalKy[0], 
					   SubsystemNbrParticles, SubsystemTotalKy, i);
				}
			      else
				{
				  sprintf (TmpEigenstateName,
					   "bosons_torus_kysym_density_n_%d_2s_%d_ky_%d_na_%d_kya_%d.%d.vec",
					   NbrParticles, KyMax, TotalKy[0], 
					   SubsystemNbrParticles, SubsystemTotalKy, i);
				}
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
	      else
		{
		  if ((EigenstateFlag == true) && (FilterKya == SubsystemTotalKy))
		    {
		      RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
						PartialDensityMatrix.GetNbrRow(), true);
		      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			TmpEigenstates[i][i] = 1.0;
		      PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
		      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
		      char* TmpEigenstateName = new char[512];
		      int MaxNbrEigenstates = NbrEigenstates;
		      if (NbrEigenstates == 0)
			MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
		      for (int i = 0; i < MaxNbrEigenstates; ++i)
			{
			  if (TmpDiag[i] > 1e-14)
			    {
			      if (Statistics == true)
				{
				  sprintf (TmpEigenstateName,
					   "fermions_torus_kysym_density_n_%d_2s_%d_ky_%d_na_%d_kya_%d.%d.vec",
					   NbrParticles, KyMax, TotalKy[0], 
					   SubsystemNbrParticles, SubsystemTotalKy, i);
				}
			      else
				{
				  sprintf (TmpEigenstateName,
					   "bosons_torus_kysym_density_n_%d_2s_%d_ky_%d_na_%d_kya_%d.%d.vec",
					   NbrParticles, KyMax, TotalKy[0], 
					   SubsystemNbrParticles, SubsystemTotalKy, i);
				}
			      TmpEigenstates[i].WriteVector(TmpEigenstateName);
			    }
			}
		      delete[] TmpEigenstateName;
		    }
		  else
		    {
		      PartialDensityMatrix.Diagonalize(TmpDiag);
		    }
		}
#else
	      PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
	      TmpDiag.SortMatrixDownOrder();
	      if (DensityMatrixFileName != 0)
		{
		  ofstream DensityMatrixFile;
		  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		  DensityMatrixFile.precision(14);
		  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKy << " " << TmpDiag[i] << endl;
		  DensityMatrixFile.close();
		}
	      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		{
		  if (TmpDiag[i] > 1e-14)
		    {
		      EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
		      DensitySum +=TmpDiag[i];
		    }
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
		    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKy << " " << TmpValue << endl;
		    DensityMatrixFile.close();
		  }		  
		if (TmpValue > 1e-14)
		  {
		    EntanglementEntropy += TmpValue * log(TmpValue);
		    DensitySum += TmpValue;
		  }
	      }
	}
      File << SubsystemNbrParticles << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << endl;
    }
  File.close();
}

