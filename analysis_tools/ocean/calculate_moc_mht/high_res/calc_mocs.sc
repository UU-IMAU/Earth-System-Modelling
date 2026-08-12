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
# Set basename of the tavg files
set basename_tavg =  mres_b.e10.B2000_CAM5.f05_t12.001

# Set directory that contains the MSF files.
set msfdir = "/projects/0/nwo2021025/archive/mres_b.e10.B2000_CAM5.f05_t12.001/ocn/msf"

# Set directory on klipdccp account where temporary empty MSF files are put when they are 
# being calculated in the queue at this moment. Then calculation is not necessary 
set msfdir_inqueue = "/home/renew1/scripts/make_moc_files/cesm/mres_b.e10.B2000_CAM5.f05_t12.001/in_msf/in_queue"

# Set directory where production run output files are located:
set tavgdir = "/projects/0/nwo2021025/archive/mres_b.e10.B2000_CAM5.f05_t12.001/ocn/hist" 

# Set directory where tool called 'calculate_msf' is located:
set tooldir = "/home/renew1/scripts/make_moc_files/cesm/code"

# You can check if the job is running by typing: mysqueue

# ----------------------------------------------------------------
# Don't change anything below 
# ----------------------------------------------------------------

# give each calculation its own id
set calcid = `awk 'FNR==1' calcid`
@ calcid_p1 = $calcid + 1
echo $calcid_p1 > calcid

if !(-e $msfdir_inqueue) mkdir -p $msfdir_inqueue

# Determine present work directory (necessary for job script)
set workdir = `pwd`

cd $tavgdir

# Determine the list of all monthly output files of the production run

ls -l ${basename_tavg}.pop.h.????-??.nc > $workdir/all_monthly_files_tmp

cd $workdir

# Make a directory for the namelist files that go into calculate_msf
#if !(-e in_msf) mkdir in_msf
#if !(-e in_msf/in_queue) mkdir -p in_msf/in_queue

awk  '{print $9}' all_monthly_files_tmp > all_monthly_files
#rm all_monthly_files_tmp

# Check the number of production run output files
set nrfiles = `cat all_monthly_files | wc -l`
echo 'nr of monthly files is: '$nrfiles

# Now loop through all files, check if there is a corresponding MSF file
# and if not make a job script that will calculate these files
  
cat $tooldir/fixed_part_job > junk
echo 'cd '$workdir >> junk

set i = 1
# This job will run on a whole node to speed up the process but because 
# it takes a lot of memory it is only possible to calculate MSF_max files 
# parrallel at a time. 

set MSF_max = 10
set MSF_id  = 1
while ($i <= $nrfiles)
  set file = `awk -v nr=$i 'FNR==nr' all_monthly_files`
  if !(-e $msfdir/MSF_${file} && $MSF_id <= $MSF_max) then
    # Check if MSF file is in $msfdir_inqueue
    if !(-e $msfdir_inqueue/MSF_${file}) then
      # determine year and month of $file  
      set datestr = `echo $file | awk 'BEGIN{FS="."}{print $4}'`
      set year = `echo $datestr | cut -b 1-4`
      set month = `echo $datestr | cut -b 5-6`
      set day = `echo $datestr | cut -b 7-8`

      # Make a namelist file in_msf_$file with help of in_template_1field  
      #sed -e 's/'fileout'/'MSF_${file}.nc'/g' $tooldir/in_template_1field_cesm > in_msf/in_msf_$file
      sed -e 's/'fileout'/'MSF_${file}'/g' $tooldir/in_template_1field_cesm > in_msf/in_msf_$file
    
      # Add input file info (=output production run) to in_msf_$file
      echo "u_file = '"$tavgdir/$file"'" >>  in_msf/in_msf_$file
      echo "v_file = '"$tavgdir/$file"'" >>  in_msf/in_msf_$file
      echo "uet_file = '"$tavgdir/$file"'" >>  in_msf/in_msf_$file
      echo "vnt_file = '"$tavgdir/$file"'" >>  in_msf/in_msf_$file
      echo '/' >>  in_msf/in_msf_$file

      # Add line to job script that will calculate the MSF file that belongs to $file
      echo $tooldir"/calculate_msf_cesm < in_msf/in_msf_"$file "&" >> junk
      sed -e 's/'calculate_moc'/'calc_moc_$year$month'/g' junk > calc_moc_all_$calcid
      echo 'Job will make MSF file for monthly file ' $file
      echo "" > $msfdir_inqueue/MSF_${file}
      @ MSF_id = $MSF_id + 1
    else
      echo 'MSF file MSF_'${file}' will be calculated by an existing job in the queue so no need to calculate it again'
    endif
  endif
  if ( ${MSF_id} == ${MSF_max} ) then
    echo "Now calculating the maximum of ${MSF_max} MSF files. After the job is done run script again to calculate the rest"
    set i = $nrfiles 
  endif
  @ i = $i + 1
end
echo 'wait' >> calc_moc_all_$calcid
echo 'mv MSF_${basename_tavg}??????.nc '$msfdir'/' >> calc_moc_all_$calcid
echo 'mv MSF_*.nc '$msfdir'/' >> calc_moc_all_$calcid

if ( ${MSF_id} == 1 ) then 
  echo 'No need to calculate any MSF file. All MSF files are in '$msfdir  
else
  echo 'Now submitting job to calculate missing  MSF files'  
  sbatch calc_moc_all_$calcid
  echo "sbatch calc_moc_all_$calcid"
  #llsubmit calc_moc_all_$calcid
  #chmod u+x calc_moc_all
  #./calc_moc_all
endif

mv MSF_*.nc $msfdir/

