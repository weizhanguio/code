export OMPI_CC=icc
CC=mpicc

##gcc ?icc?  -ansi_alias -prefetch
app:  AVX_mpi_3D.c
	$(CC)    AVX_mpi_3D.c  -o app  -O2  -xavx 

clean:
	rm -rf app  *.o    *.*~  *.e*  *.txt  *~
