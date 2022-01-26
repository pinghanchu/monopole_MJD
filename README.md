# monopole_MJD

1.Modify from the $GEANT4DIR/examples/extended/exoticphysics/monopole/.

2.Output the stepping parameters.

3.Use Majorana Demonstrator detectors.

4.mkdir ./monopole-build/

cd monopole-build/

cmake -DGeant4_DIR=/mjsw/GEANT/install/lib64/Geant4-10.4.3/ ~/work/g4work/monopole_MJD/

make -j4

#Analysis
1. ./monopole monopole.in
This generates monopole.root

2. Need to modify the monopole mass in ./src/G4MonopolePhysics.cc and re-compile before generating Monte Carlo
Modify the line "fMonopoleMass = 1e12*GeV;"

3. ./analysis/run.pl 1e12
1e12 is the monopole mass used in G4MonopolePhysics.cc
The script will run replace.pl and generate Monte Carlo with the energy from 1e-5 GeV to 1e10 GeV.

4. ./analysis/replace.pl $erg
$erg is the monopole energy. The script will modify the monopole energy in monopole.in.
./analysis/replacemass.pl $mass
$mass is the monopole mass.

5. ./analysis/build.pl 1e12GeV
1e12GeV is the monopole mass.
The script will generate data.csv including information of each step.

6. ./analysis/fill.cc 
The macro will generate csv data from monopole.root.

7. ./analysis/gatify.py
The script will re-construct steps as events. ./data/event_[mass]_[energy].h5/csv
Need to run python in the default environment of CORI and ">module load python".

#Monte Carlo
The unit of Monte Carlo: time(ns), energy(MeV), length(mm)

data1: vertical monopole beam
data2: hemispherical monopole beam