# Algorithm Input and Output Data Definition (IODD)

# Input data

The input data for the multi-parameter retrieval includes all channels at
horizontal and vertical polarization and their uncertainties coming from the
L1b processor. In addition, time and location is required for calculation of
field of resampled TBs. The evaluation of the algorithm is planned on the
C-band footprints, with all other channels being resampled to it. The
uncertainties are assumed Gaussian in the retrieval. 

## Output data

The output data consist of the parameters, namely {term}`WSP`, {term}`TWV`,
{term}`CLW`, {term}`SST`, {term}`IST`, {term}`SIC`, {term}`MYIF` and
{term}`SIT` as well as time and location for the resampling. The uncertainties are assumed
Gaussian in the retrieval and are provided along with the parameters. 

## Auxiliary data
Auxillary data includes the {term}`TEC` for the transformation to account for
the polarization rotation as part of the atmospheric model. {term}`ECMWF`
surface analysis data is used as background values for the retrieval. The
variables used are {term}`WSP`, {term}`TWV`, {term}`CLW`, {term}`T2M`,
{term}`TSK`. 


## Ancillary data
No ancillary data is used in the retrieval.
