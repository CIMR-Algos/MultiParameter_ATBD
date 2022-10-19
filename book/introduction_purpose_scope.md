# Introduction, purpose and scope

This is the {term}`ATBD` of the multi parameter retrieval. The method is based
on the works of {cite}`Pedersen1991,Scarlat2017,Scarlat2018,Scarlat2020` and
will be described in this document. It is applied to {term}`CIMR` L1b data which
is resampled to a common footprint for all frequencies (we call L1R). It is
using all {term}`CIMR` frequency channels, namely 1.4, 6.9, 10.7, 18.7, and 36.5
GHz and {term}`ECMWF` Analysis as input. The output of this multi parameter
retrieval is in the same resolution format, i.e. L2R. It is physically
consistent and can be used in turn as a priori for the other retrieval
parameters. 

This document is describing the algorithm and processing steps of the L2R multi
parameter retrieval product. The document is intended for the {term}`CIMR` users and
interested parties. It is not intended to replace a product user guide. The
algorithm is implemented in the multi parameter retrieval software package which
is distributed in the "algorithm" directory of this repository.
