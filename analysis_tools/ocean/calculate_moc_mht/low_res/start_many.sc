#!/bin/csh -f 

set N = $argv[1]

echo 'nr of times N is: ' $N

set i=1
while ($i <= $N)
   echo ""
   echo "$i th time of $N"
   echo ""
  ./calc_mocs.sc
  @ i = $i + 1
end


