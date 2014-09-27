export OMPI_CC=icc
CC=mpicc
##gcc ?icc?  -ansi_alias -prefetch
app:  AVX_mpi_1D.c
	$(CC)    AVX_mpi_1D.c  -o app  -O2 -vec-report0  -xavx

clean:
	rm -rf app  *.o    *.*~  *.e*  *.txt  *~
