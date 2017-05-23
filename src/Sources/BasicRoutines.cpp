/**
 * @file        BasicRoutines.cpp
 * @author      Jiri Jaros and Vojtech Nikl\n
 *              Faculty of Information Technology \n
 *              Brno University of Technology \n
 *              jarosjir@fit.vutbr.cz
 *
 * @brief       The implementation file with basic routines
 *
 * @version     2015
 * @date        02 March 2015, 14:09 (created) \n
 *              02 March 2015, 16:22 (revised)
 *
 * @detail
 * This implementation file with basic routines and simulation parameters
 */

#include <getopt.h>
#include <string>
#include <cmath>
#include <iostream>

#include "BasicRoutines.h"
#include "MaterialProperties.h"



/**
 * Print parameters
 */
void TParameters::PrintParameters() const
{
  if (batchMode)
  {
    /*
    printf("%d;%ld;%ld;%ld;%2.6f;%s;",worldSize,edgeSize, nIterations,
                                    diskWriteIntensity, airFlowRate, 
                                    materialFileName.c_str());
                                    */
    cout << "MPIProcs:" << worldSize << endl;
    cout << "Domain:" << edgeSize << endl;
    cout << "Iterations:" << nIterations << endl;
    cout << "DiskWrite:" << diskWriteIntensity << endl;
    cout << "AirFlow:" << airFlowRate << endl;
    cout << "InputFile:" << materialFileName.c_str() << endl;
    cout << "Imbalance threshold: " << std::fixed << threshold << endl;

  }
  else
  {
    printf(".......... Parameters of the simulation ...........\n");
    printf("Processes           : %d \n", worldSize);
    printf("Domain size:          %ldx%ld\n", edgeSize, edgeSize);
    printf("Number of iterations: %ld\n",     nIterations);
    printf("Disk write intensity: %ld\n",     diskWriteIntensity);
    printf("Air flow rate       : %2.6f\n",   airFlowRate);
    printf("Input file name     : %s \n",     materialFileName.c_str());
    printf("Output file name    : %s \n",     outputFileName.c_str());
    printf("Mode                : %d \n",     mode);
    printf("Object size         : %d \n",     objDim);
    printf("...................................................\n\n");
  }
} // end of PrintParameters
//------------------------------------------------------------------------------


/**
 * Parsing arguments from command line
 * @param [in] argc
 * @param [in] argv
 * @param [out] parameters - returns filled struct
 */
void ParseCommandline(int argc, char *argv[], TParameters &parameters)
{
  int c;

  bool n_flag = false;
  bool w_flag = false;
  bool i_flag = false;
  bool m_flag = false;
  bool T_flag = false;
  bool M_flag = false;

  string temp, xs,ys;

  while ((c = getopt (argc, argv, "n:w:a:dvi:o:bm:ps:t:XT:M:")) != -1)
  {
    switch (c)
    {
      case 'n':
        parameters.nIterations = atol(optarg);
	      n_flag = true;
        break;

      case 'w':
        parameters.diskWriteIntensity = atol(optarg);
	      w_flag = true;
        break;

      case 'a':
        parameters.airFlowRate = atof(optarg);
        break;

      case 'd':
        parameters.debugFlag = true;
        break;

      case 'm':
	      parameters.mode = atol(optarg);
        m_flag = true;
        break;

      case 'v':
	      parameters.verificationFlag = true;
        break;

      case 'i':
        i_flag = true;
        parameters.materialFileName.assign(optarg);
        break;

      case 'o':
        parameters.outputFileName.assign(optarg);
        break;
      
      case 'b':
        parameters.batchMode = true;
        break;

      case 'p':
        parameters.useParallelIO = true;
        break;

      case 's':

        temp.assign(optarg);
        // xs = temp.substr(0, temp.find(':'));
        parameters.objDim = stoi(temp);

        // ys = temp.substr(temp.find(':')+1, temp.size());
        // parameters.objDimY = stoi(ys);
        break;

      case 't':

        parameters.balancePeriod = atol(optarg);
        break;

      case 'X':
        parameters.balance = true;
        break;

      case 'T':
        parameters.threshold = atof(optarg);
        T_flag = true;
        break;

      case 'M':
        parameters.multiply = atof(optarg);
        M_flag = true;
        break;

      default:
        fprintf(stderr,"Wrong parameter!\n");
        PrintUsageAndExit();
    }
  } // while

  MPI_Comm_size(MPI_COMM_WORLD, &(parameters.worldSize));
  // output enable flag
  parameters.ioEnabled = (parameters.outputFileName != "");
  if(! T_flag)
    parameters.threshold = 1.5;

  if(! M_flag)
    parameters.multiply = 1;


  if (!(n_flag && i_flag && w_flag && m_flag) || 
      !(parameters.mode >= 0 && parameters.mode <= 2))
  {
    PrintUsageAndExit();
  }
} // end of ParseCommandline
//------------------------------------------------------------------------------


/**
 * Verify the difference between two corresponding gridpoints in the sequential
 * and parallel version. Print possible problems
 * @param [in] tempsSeq
 * @param [in] tempsPar
 * @param [in] parameters
 * @param [in] epsilon
 * @return true if no problem found
 */
bool VerifyResults(const float      *seqResult,
                   const float      *parResult,
                   const TParameters parameters,
                   const float       epsilon)
{
  for (size_t i = 0; i < parameters.edgeSize * parameters.edgeSize; i++)
  {
    if (fabs(parResult[i] - seqResult[i]) > epsilon)
    {
      printf("Error found at position -> difference: ");
      printf("[%ld, %ld] -> %e\n", i / parameters.edgeSize, i % parameters.edgeSize,
                                    parResult[i] - seqResult[i]);
      return false;
    }
  } // for

  return true;
} //end of VerifyResults
//------------------------------------------------------------------------------


/**
 * Print usage and exit.
 */
void PrintUsageAndExit()
{
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"Mandatory arguments:\n");
  fprintf(stderr,"  -m [0-2]    mode 0 - run sequential version\n"); 
  fprintf(stderr,"              mode 1 - run parallel version\n");
  fprintf(stderr,"  -n number of iterations\n");
  fprintf(stderr,"  -w disk write intensity (how often)\n");
  fprintf(stderr,"  -i material hdf5 file\n");

  fprintf(stderr,"Dynamic Load Balancing - optional\n\n");
  fprintf(stderr,"  -X set balancing ON (default OFF)\n");
  fprintf(stderr,"  -t number of iterations in one balancing period \n");
  fprintf(stderr,"     imbalance detection will run once per t iterations (default nIterations / 10) \n" );
  fprintf(stderr,"  -s object size (parallle mode), format x:y (default 8:8)\n\n" );

  fprintf(stderr,"Optional arguments:\n");
  fprintf(stderr,"  -o output hdf5 file\n");
  fprintf(stderr,"  -a air flow rate (values within <0.5, 0.0001> make sense\n");
  fprintf(stderr,"  -d set debug mode (compare results from seq and par version and write them to cout)\n");
  fprintf(stderr,"  -v verification mode (compare results of seq and par version)\n");
  fprintf(stderr,"  -p parallel I/O mode\n");
  fprintf(stderr,"  -b batch mode - output data in CSV format\n");
  fprintf(stderr,"  -M delay multiplier - float\n");
  fprintf(stderr,"  -T balancing threshold - float\n");

  
  exit(EXIT_FAILURE);
} // end of print usage
//------------------------------------------------------------------------------


/**
 *  Print array content
 * @param [in] data
 * @param [im] size
 */
void PrintArray(const float *data, const size_t edgeSize)
{
  for (size_t i = 0; i < edgeSize; i++)
  {
    printf("[Row %ld]: ", i);
    for (size_t j = 0; j < edgeSize; j++)
    {
      printf("%e, ", data[i * edgeSize + j]);
    }
    printf("\n");
  }
} // end of PrintArray
//------------------------------------------------------------------------------
