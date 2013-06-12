#!/bin/bash

touch all_trajectories.txt

para=1
while [ $para -le 128 ]; do
    cd $para
      
#      pwd
#      clear all;
      cp trajectories.txt ../
      cd ..
      cat trajectories.txt >> all_trajectories.txt
      rm trajectories.txt
   let para=para+1
   done
