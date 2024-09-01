# ENUBET

## Contains 

FIX_MARCHESINI.C pass from a Run(numero)_list.root file to a completly fixed (e, z, phi, LG, HG, timestamp) thank to the mapping file (here called mappingMarchesini.txt), the MappingScript.h and the calibration given by the test beam (calibrationMarchesini)
  The only sector to modify by defining the parameters
  - p is the file generated as mean and std of the pedestal
  - filename is the input file given by Janus
  - output is the name of the final output file
  - fileInput is the calibration of channels given by board/channel/weight
  - numberOfBoards is the number of boards
  - thresholdTrigger is the time coincidence in mu s to impose
  - sigma is the number of std of the pedestal to which it eliminates the value
  - numberOfAnodes is the number of channels connected to the board
  - majority is a possible secondary majority with respect to the one of Janus

The code works by first doing a gaussian fit of all the data, only events above a certain value mean+NumberOfSigma*sigmaFit, then by mapping and equalizing the detector, at the end we need to pass from a board based file to an event based and that is done by chosing the time coincidence window 
- ANGLE.C from a r, z, phi, LG and HG file execute a linear 3d fit returning a root file with azimuth, elevation, errror on both, direction X, direction Y, direction Z and relative errors, NDF and the MinFCN found
- sum0.C execute the sum over each single channel and return a .txt file with the hit map values
- sum1.C does a plot with the output file of sum0.C but in general with every .txt file composed similarly, could be also the mean value, sigma of the fit...
- fulldistr.C generate a .root file with events of cosmic muons based on a momentum dependent probability distribution function (check the report for more details)
- merge.C can be used to merge multiple files given by the output of the FIX_MARCHESINI.C program to have higher statictics
- fG.C using the mappingMarchesiniActive.txt file reorganize the templateAnlysis.root file output of the Geant4 simulation in a (r, z, phi, energy) file
