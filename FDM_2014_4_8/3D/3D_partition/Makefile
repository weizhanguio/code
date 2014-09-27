export OMPI_CC=icc
CC=mpicc

##gcc ?icc?  -ansi_alias -prefetch
app:  mix_3D_partition.c
	$(CC)    mix_3D_partition.c  -o app  -O2 -openmp -xavx

clean:
	rm -rf app  *.o    *.*~  *.e*  *.txt  *~
