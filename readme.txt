mkdir B1-build
cd B1-build/
cmake -DGeant4_DIR=/mjsw/GEANT/install/lib64/Geant4-10.4.3/ ~/work/g4work/B1/
make -j4
./exampleB1 
