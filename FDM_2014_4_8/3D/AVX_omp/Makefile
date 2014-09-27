CC=icc

app:  AVX_omp.c
	$(CC)    AVX_omp.c  -o app  -O2 -openmp   -xavx

clean:
	rm -rf app  *.o    *.*~  

