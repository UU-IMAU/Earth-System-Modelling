# Calculate MOC and MHT from POP output files

This directory contains scripts and code that calculates the Meridional Overturning Circulation (MOC)
in Sverdrup and the mean Meridional Heat Transport (MHT) in PetaWatts both for the Atlantic,
the Indo Pacific and the World, from fields in the monthly mean outputfile 
of the Parallel Ocean Program (POP) model. 

The code needs masks of the Atlantic and IndoPacific as well as the horizontal grid and depth levels
that is why there is a distinction between low (1 degree) and high (0.1 degree) resolution POP output.
