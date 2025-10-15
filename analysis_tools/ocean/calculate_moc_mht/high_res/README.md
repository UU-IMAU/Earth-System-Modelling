# Introduction

`calc_moc.sc` is a **C-shell script** that calculates an *MSF (Meridional Stream Function)* file from a **monthly mean POP output file**.

The generated MSF file contains the following fields:

| Field  | Description | Units |
|---------|--------------|--------|
| **MHTG**  | Mean meridional heat transport (Global, i.e., total ocean) | PetaWatts |
| **MHTA**  | Mean meridional heat transport (Atlantic) | PetaWatts |
| **MHTIP** | Mean meridional heat transport (Indo-Pacific) | PetaWatts |
| **TMTG**  | Meridional overturning streamfunction (Global, i.e., total ocean) | Sverdrups (Sv) |
| **TMTA**  | Meridional overturning streamfunction (Atlantic) | Sverdrups (Sv) |
| **TMTIP** | Meridional overturning streamfunction (Indo-Pacific) | Sverdrups (Sv) |

---

# Usage

`calc_moc.sc` calls the Fortran program [`calculate_msf_cesm.f`](code/calculate_msf_cesm.f), which performs the actual MOC (Meridional Overturning Circulation) calculation.

Before using `calc_moc.sc`, you must **compile** the Fortran program.

### 1. Load required modules

```
module purge
module load 2023
module load intel/2023a
module load netCDF/4.9.2-iimpi-2023a
module load netCDF-Fortran/4.6.1-iimpi-2023a
```

### 2. Compile the program

```
cd code
ifort -O3 -extend-source -convert big_endian \
  -o calculate_msf_cesm calculate_msf_cesm.f \
  -lnetcdf -lnetcdff
```

### 3. Run the script

After compilation, you can easily run the tool from the command line, for example:

```
./calc_moc.sc /projects/0/prace_imau/prace_2013081679/cesm1_0_4/f05_t12/spinup_pd_maxcores_f05_t12/OUTPUT/ocn/hist/monthly/spinup_pd_maxcores_f05_t12.pop.h.0300-01.nc
```

The tool calculates the MSF fields based on the UVEL and VVEL fields in the specified POP output file.

### 4. Monitoring the Job

You can check if the job is running using:

```
mysqueue
```

To monitor progress, follow the job’s output file (e.g., slurm-123456.out):

```
tail -f slurm-123456.out
```

### 5. Running Multiple Jobs

If you want to run calc_moc.sc multiple times in sequence and calculate all monthly files of specific years you can use:

```
sbatch calc_moc_many.sc
```

