#!/bin/csh -f
# ----------------------------------------------------------------
#
# This script checks which MSF files are not in directory $msfdir (defined below) yet
# and calculates them via a job script. It can calculate $MSF_max  (defined below) files 
# at the same time. The new MSF files will automatically be put in $msfdir.
#
# An MSF file is a file containing MOC and Heat Transport
# information based on a monthly mean outputfile of a POP run.
# More information about the MSF files can be found in the README file.
#
# author: Michael Kliphuis
#
# ----------------------------------------------------------------
# User section 
# ----------------------------------------------------------------
#
# Set basename of the monthly mean POP outputfiles, the so called tavg files
set basename_tavg = t.x1_SAMOC_flux.

# Set directory where the calculated MSF files will be put.
set msfdir = "msf"   # this can also be a directory like: "/projects/0/uus20475/pop/gx1v6/x1_SAMOC_flux/tavg"

# Set directory on Snellius account where temporary and empty MSF files will be put as soon as a job is in the queue.
# As long as a job is in the queue then calculating this MSF file again is not needed.
set msfdir_inqueue = "in_msf/in_queue"

# Set directory where the monthly mean POP outputfiles are located:
set tavgdir = "tavg"

# Set directory where tool called 'calculate_msf_1deg' and other code files are located:
set tooldir = "code"

# Determine present work directory. This script automatically generates a job script that will be put here.
set workdir = `pwd`

# This job will run on a whole node to speed up the process but because 
# it takes memory it is better to calculate $MSF_max files 
# parrallel at a time instead of 192 (one per processor) on a Genoa node. 
set MSF_max = 100   # should be <= 192 i.e. number of cores on a Genoa node

# ----------------------------------------------------------------
# Don't change anything below, only if really needed 
# ----------------------------------------------------------------

# create necessary directories
if !(-e $msfdir) mkdir -p $msfdir
if !(-e $msfdir_inqueue) mkdir -p $msfdir_inqueue

# Determine present work directory (necessary for job script)
cd $tavgdir

# Determine the list of all monthly output files of the production run

ls -l ${basename_tavg}*[^hdr][^nc] > $workdir/all_monthly_files_tmp

cd $workdir

# give each calculation its own id
set calcid = `awk 'FNR==1' calcid`
@ calcid_p1 = $calcid + 1
echo $calcid_p1 > calcid

awk  '{print $9}' all_monthly_files_tmp > all_monthly_files
rm all_monthly_files_tmp

# Check the number of production run output files
set nrfiles = `cat all_monthly_files | wc -l`
echo 'nr of monthly files is: '$nrfiles

# Now loop through all files, check if there is a corresponding MSF file
# and if not make a job script that will calculate these files
  
cat $tooldir/fixed_part_job_1deg > junk
echo 'cd '$workdir >> junk

set i = 1

# This job will run on a whole node to speed up the process but because 
# it takes memory it is better to calculate $MSF_max files 
# parrallel at a time instead of 192 (one per processor) on a Genoa node. 
# If more than $MSF_max files need to be calculated then wait until these
# $MSF_max files are made and then run this script again to calculate the next ones.

set MSF_id = 1
while ($i <= $nrfiles)
  set file = `awk -v nr=$i 'FNR==nr' all_monthly_files`
  if !(-e $msfdir/MSF_${file}.nc && $MSF_id <= $MSF_max) then
    # Check if MSF file is in $msfdir_inqueue, in that case it is not needed to calculate it again
    if !(-e $msfdir_inqueue/MSF_${file}.nc) then
      # determine year and month of $file  
      set datestr = `echo $file | awk 'BEGIN{FS="."}{print $4}'`
      set year = `echo $datestr | cut -b 1-4`
      set month = `echo $datestr | cut -b 5-6`
      set day = `echo $datestr | cut -b 7-8`

      # Make a namelist file called in_msf_$file with help of in_template_1field_1deg  
      # Replace strings 'fileout' and 'tooldir' respectively by strings MSF_${file}.nc and $tooldir
      sed -e 's/'fileout'/'MSF_${file}.nc'/g' -e 's/'tooldir'/'$tooldir'/g' $tooldir/in_template_1field_1deg > in_msf/in_msf_$file
    
      # Add input file info to in_msf_$file
      echo "u_file = '"$tavgdir/$file"'"   >>  in_msf/in_msf_$file
      echo "v_file = '"$tavgdir/$file"'"   >>  in_msf/in_msf_$file
      echo "uet_file = '"$tavgdir/$file"'" >>  in_msf/in_msf_$file
      echo "vnt_file = '"$tavgdir/$file"'" >>  in_msf/in_msf_$file
      echo '/' >>  in_msf/in_msf_$file

      # Add a line to the job script that will calculate the MSF file that belongs to $file
      echo $tooldir"/calculate_msf_1deg < in_msf/in_msf_"$file "&" >> junk
      
      sed -e 's/'calculate_moc'/'calc_moc_$year$month'/g' junk > calc_moc_all_$calcid
      echo 'Job will make MSF file for monthly file ' $file
      echo "" > $msfdir_inqueue/MSF_${file}.nc
      @ MSF_id = ${MSF_id} + 1
    else
      echo 'MSF file MSF_'${file}'.nc will be calculated by an existing job in the queue so no need to calculate it again'
    endif
  endif

  if ( ${MSF_id} > ${MSF_max} ) then
    echo "Now calculating the maximum of ${MSF_max} MSF files. After the job is done run script again to calculate the rest"
    break
  endif
  @ i = $i + 1
end

echo 'wait' >> calc_moc_all_$calcid
echo 'mv MSF_*.nc '$msfdir'/' >> calc_moc_all_$calcid
echo 'rm junk all_monthly_files' >> calc_moc_all_$calcid

echo "Moved all created MSF files to directory $msfdir"

if ( ${MSF_id} == 1 ) then 
  echo 'No need to calculate any MSF file. All MSF files are in '$msfdir  
else
  echo 'Now submitting job to calculate missing  MSF files'  
  sbatch calc_moc_all_$calcid
  # uncomment lines below if you want to calculate them interactively
  #chmod u+x calc_moc_all_$calcid
  #./calc_moc_all_$calcid
endif



