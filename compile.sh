clear
cd ~/2DFit_simulations																																													
rm -f 2DFit64																																																									# delete old binary
cd ~/2DFit_simulations/obj
rm -f 2DFit64																																																									# delete old binary
cd ~/2DFit_simulations
make -j30                                                                                                                       #compile
cp ~/2DFit_simulations/obj/2DFit64  ~/2DFit_simulations                                      #copy exectutable
cd ~/2DFit_simulations

for i in {1..1}; do
mpiexec -n 72 ~/2DFit_simulations/2DFit64  
done
																																																																																# test mpi

