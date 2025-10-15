# Introduction 

  `calc_moc.sc` is a csh script that calculates a so called MSF file for a monthly mean POP outputfile.

  The generated MSF file contains the following fields:

    MHTG    Mean meridional heat transport Global i.e. world in PetaWatts
    MHTA    Mean meridional heat transport Atlantic in PetaWatts
    MHTIP   Mean meridional heat transport IndoPacific in PetaWatts
    TMTG    Meridional overturning streamfunction Global i.e. world in Sv
    TMTA    Meridional overturning streamfunction Atlantic in Sv
    TMTIP   Meridional overturning streamfunction IndoPacific in Sv

# Usage

  `calc_moc.sc` calls the fortran program `calculate_msf_cesm.f` which does the actual calculation of the MOC.
  Before using `calc_moc.sc` you first need to compile this fortran program. Do this by first loading the needed modules by typing:

module purge
module load 2023
module load intel/2023a
module load netCDF/4.9.2-iimpi-2023a
module load netCDF-Fortran/4.6.1-iimpi-2023a

Then do the actual compilation by typing:

`cd code`
`ifort -O3 -extend-source -convert big_endian -o calculate_msf_cesm calculate_msf_cesm.f -lnetcdf -lnetcdff`

  After that usage is easy, just type on the commandline e.g.:

  `./calc_moc.sc /projects/0/prace_imau/prace_2013081679/cesm1_0_4/f05_t12/spinup_pd_maxcores_f05_t12/OUTPUT/ocn/hist/monthly/spinup_pd_maxcores_f05_t12.pop.h.0300-01.nc`

  The tool then calculates the MSF fields based on the UVEL and VVEL fields in the POP output file `spinup_pd_maxcores_f05_t12.pop.h.0300-01.nc`
  
  Check if the job is running by typing:  
  
  `mysqueue`

  You can follow the evolution of the job by 'tailing' the `slurm-xxxxxx.out` (e.g. `slurm-123456.out`)
  job log file. Do this by typing:  
  
  `tail -f slurm-123456.out`

  If you want to run `./calc_moc.sc` many times after each other then type e.g.

  `./calc_moc_many.sc`

