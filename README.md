# monopole_MJD

1.Modify from the $GEANT4DIR/examples/extended/exoticphysics/monopole/.

2.Output the stepping parameters.

3.Use Majorana Demonstrator detectors.

4.mkdir ./monopole-build/
cd monopole-build/
cmake -DGeant4_DIR=/mjsw/GEANT/install/lib64/Geant4-10.4.3/ ~/work/g4work/monopole_MJD/
make -j4
./monopole
