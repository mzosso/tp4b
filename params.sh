#!/bin/bash


# set the include path to find TRandom.h


# compile the tp4b.cc file
#g++ -std=c++14 -I$ROOTSYS/include tp4b.cc -o tp4b.exe `root-config --glibs`
# Define arrays for n_d, T, and V values
#n_ds=(50 100 200 300 350)
#n_ds=(50 100)
#Ts=(100 300)
#Vs=(300 400)
Ts=(110 200 430 500 550 300)
Vs=(300 1000 1500 2000 600)
#Areas=(10000 1000000)
Areas=(0.0000035 0.0000065 0.0000045 0.0000075 0.0000055)
Capacities=(0.00000000004 0.00000000003 0.00000000002 0.00000000005 0.00000000001)
Noises=(1 2 2.5 3 3.2)
Thresholds=(2 3 3.5 4 5 7 8)
Pourcentages=(0.2 0.3 0.4 0.5 0.6 0.7 0.8)

#200 microm, 600 v, 300 K, area 5.5 mm^2, amplitude 6 mV, sigma= 3.5 mV  

#capacidad, ruido, threshold 
# Loop over values for n_d
#for n_d in "${n_ds[@]}"
#do
  
 
  # Update n_d, T, and V values in config1.txt
     # sed -i '' -e "2s/.*/$n_d/" config1.txt
      
     # root -q tp4b.cc >> "output_n_d.txt" 
      #./tp4b.exe >> output_n_d_
      # Run tp4b.exe and save the output to a file
      
#done

#line 14: 0 if not CFD, 1 if CFD, line 16: 0 if capacity from area and thickness, 1 if not
sed -i '' -e '14s/.*/0/' config1.txt
sed -i '' -e '16s/.*/0/' config1.txt
for area in "${Areas[@]}"
do
  
 
  # Update n_d, T, and V values in config1.txt
      sed -i '' -e "1s/.*/$area/" config1.txt
      echo $area
      
      root -q tp4b.cc >> "output_area.txt" 
      #./tp4b.exe >> output_n_d_
      # Run tp4b.exe and save the output to a file
      
done

 # Loop over values for T
   for T in "${Ts[@]}"
   do
      sed -i '' -e "7s/.*/$T/" config1.txt
      root -q tp4b.cc >>output_T.txt 2>&1
   done

    # Loop over values for V
   for V in "${Vs[@]}"
     do
        sed -i '' -e "4s/.*/$V/" config1.txt
 
        root -q tp4b.cc >>output_V.txt 2>&1
      
       # Extract the time resolution from the output file and print it to the terminal
   
     done
     sed -i '' -e '16s/.*/1/' config1.txt
     for C in "${Capacities[@]}"
     do
        sed -i '' -e "15s/.*/$C/" config1.txt
 
        root -q tp4b.cc >>output_C.txt 2>&1
      
       # Extract the time resolution from the output file and print it to the terminal
   
     done
  sed -i '' -e '16s/.*/0/' config1.txt
   for noise in "${Noises[@]}"
     do
        sed -i '' -e "13s/.*/$noise/" config1.txt
 
        root -q tp4b.cc >>output_noise.txt 2>&1
      
       # Extract the time resolution from the output file and print it to the terminal
   
     done


sed -i '' -e '14s/.*/0/' config1.txt
    for threshold in "${Thresholds[@]}"
     do
        sed -i '' -e "12s/.*/$threshold/" config1.txt
        root -q tp4b.cc >>output_threshold.txt 2>&1
      
       # Extract the time resolution from the output file and print it to the terminal
   
      done

  sed -i '' -e '14s/.*/1/' config1.txt
    for pourcentage in "${Pourcentages[@]}"
     do
        sed -i '' -e "17s/.*/$pourcentage/" config1.txt
        root -q tp4b.cc >>output_porcentage.txt 2>&1
      
       # Extract the time resolution from the output file and print it to the terminal
   
      done
  


