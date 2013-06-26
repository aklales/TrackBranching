#!/bin/bash

# The following script will 
# (1) create 'x' number of directories, equal to number of cores #      we want to use
# (2) each directory will contain a subset of the particles
# (3) and also a bsub file
# (4) finally, it will be submitted to the queue

para=1
while [ $para -le 4 ]; do
   mkdir $para
   cp Potential.mat ./$para
   cp betterbinz.m ./$para
   cp Propagate_interpolate.m ./$para
   cd $para
      pwd

      clear all;
      echo "#!/bin/sh" > gobranch.bsub
      echo "#BSUB -n 1" >> gobranch.bsub
      echo "#BSUB -J branchjob${para}" >> gobranch.bsub
      echo "#BSUB -o branch_lsf.out" >> gobranch.bsub
      echo "#BSUB -e branch_lsf.err" >> gobranch.bsub
      echo "#BSUB -q short_serial" >> gobranch.bsub
      echo "module load math/matlab-R2012b" >> gobranch.bsub
      echo 'matlab -nosplash -nodesktop -r "Propagate_interpolate"' >> gobranch.bsub

#make -f labmake.make

      bsub < gobranch.bsub

  cd ..
  let para=para+1
  done
