#!/bin/bash


# set the include path to find TRandom.h


# compile the tp4b.cc file
#g++ -std=c++14 -I$ROOTSYS/include tp4b.cc -o tp4b.exe `root-config --glibs`
# Define arrays for n_d, T, and V values
#n_ds=(50 100 200 300 350)
n_ds=(50 100)
Ts=(6 45 110 200 300 430)
Vs=(300 400 500 550 600)

# Loop over values for n_d
for n_d in "${n_ds[@]}"
do
  
 
  # Update n_d, T, and V values in config1.txt
      sed -i '' -e "2s/.*/$n_d/" config1.txt
      
      root -q tp4b.cc >> "output_n_d.txt" 2>&1
      #./tp4b.exe >> output_n_d_
      # Run tp4b.exe and save the output to a file
      
done

# # Loop over values for T
#   for T in "${Ts[@]}"
#   do
#       sed -i '' -e "7s/.*/$T/" config1.txt
#     root -q tp4b.cc >>output_T.txt 2>&1
#   done

#    # Loop over values for V
#   for V in "${Vs[@]}"
#     do
#       sed -i '' -e "4s/.*/$V/" config1.txt
 
#      root -q tp4b.cc >>output_V.txt 2>&1
      
#       # Extract the time resolution from the output file and print it to the terminal
   
#     done
  


