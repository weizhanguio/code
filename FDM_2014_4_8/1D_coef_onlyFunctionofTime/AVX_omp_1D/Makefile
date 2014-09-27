CC=icc

app:  AVX_omp_1D.c
	$(CC)    AVX_omp_1D.c  -o app  -O2 -vec-report0 -openmp -xavx

clean:
	rm -rf app  *.o    *.*~  *.e*  *.txt 
