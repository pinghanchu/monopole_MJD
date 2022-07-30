# monopole_MJD

1. This package was modify from the $GEANT4DIR/examples/extended/exoticphysics/monopole/.

2. Output the stepping parameters.

3. Use Majorana Demonstrator detectors.

# Compile monopole_MJD
1. $mkdir ./monopole-build/

2. $cd ./monopole-build/

3. $cmake -DGeant4_DIR=/mjsw/GEANT/install/lib64/Geant4-10.4.3/ ../

4. $make

# Excute monopole_MJD

1. Need copy ./Detectorposition.txt to ./monopole-build/

2. $./monopole monopole.in
This generates monopole.root

3. Need to modify the monopole mass in ./src/G4MonopolePhysics.cc and re-compile before generating Monte Carlo
Modify the line "fMonopoleMass = 1e12*GeV;"

## Script
1. $./analysis/run.pl 1e12
1e12 is the monopole mass used in G4MonopolePhysics.cc
The script will run replace.pl and generate Monte Carlo with the energy from 1e-5 GeV to 1e10 GeV.

2. $./analysis/replace.pl $erg
$erg is the monopole energy. The script will modify the monopole energy in monopole.in.
./analysis/replacemass.pl $mass
$mass is the monopole mass.

# Analysis
1. $./analysis/analysis.cc 
   $./analysis/analysis.pl
This script reconstructed and analyzed the raw data (monopole.root)

# Monte Carlo
The unit of Monte Carlo: time(ns), energy(MeV), length(mm)