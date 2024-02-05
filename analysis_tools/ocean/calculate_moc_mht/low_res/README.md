# Introduction 

  `calc_mocs.sc` is a csh script that calculates the so called MSF files for all monthly mean POP outputfiles in
  a certain directory `tavgdir` (defined in `calc_mocs.sc`) containing monthly mean binary POP outputfiles.

  The generated MSF files contain the following fields:

    MHTG    Mean meridional heat transport Global i.e. world in PetaWatts
    MHTA    Mean meridional heat transport Atlantic in PetaWatts
    MHTIP   Mean meridional heat transport IndoPacific in PetaWatts
    TMTG    Meridional overturning streamfunction Global i.e. world in Sv
    TMTA    Meridional overturning streamfunction Atlantic in Sv
    TMTIP   Meridional overturning streamfunction IndoPacific in Sv

# Usage

  `calc_mocs.sc` calls the fortran program `calculate_msf_1deg.f` which does the actual calculation of the MOC.
  Before using `calc_mocs.sc` you first need to compile this fortran program. Do this by typing:

  `cd code`

  Then on the national supercomputer Snellius do the compilation by typing: 
  
  `gfortran -O3 -fconvert=big-endian -I/sw/arch/RHEL8/EB_production/2022/software/netCDF-Fortran/4.6.0-gompi-2022a/include -o calculate_msf_1deg calculate_msf_1deg.f -lnetcdf -lnetcdff`

  After that usage is easy, just type on the commandline:

  `./calc_mocs.sc`

  The tool then makes a list of all the monthly mean binary POP outputfiles e.g. `t.x1_SAMOC_flux.yyyymm`
  in the directory `tavgdir`
  
  For each file in this list it checks if there is a corresponding MSF file and then makes it if it is not
  there. More precise: it automatically generates a job script that calculates up to `MSF_max` (defined in `calc_mocs.sc`) 
  files at the same time.

  Check if the job is running by typing:  
  
  `mysqueue`

  You can follow the evolution of the job by 'tailing' the `slurm-xxxxxx.out` (e.g. `slurm-123456.out`)
  job log file. Do this by typing:  
  
  `tail -f slurm-123456.out`

  If you want to run `./calc_mocs` many times after each other then type e.g.

  `./start_many.sc 50`

  this will start it 50 times.
