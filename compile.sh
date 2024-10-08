clear
cd ~/2DFit																																													
rm -f 2DFit64																																																									# delete old binary
cd ~/2DFit/obj
rm -f 2DFit64																																																									# delete old binary
cd ~/2DFit
make -j30                                                                                                                       #compile
cp ~/2DFit/obj/2DFit64  ~/2DFit                                     #copy exectutable
cd ~/2DFit

for i in {1..1}; do
mpiexec -n 72 ~/2DFit/2DFit64  
done
																																																																																# test mpi

