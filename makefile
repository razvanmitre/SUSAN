all: serial mpi omp 

serial: susan.c
	gcc susan.c -lm -o susan_serial

mpi: susan_mpi.c
		

omp: susan_omp.c
	gcc -fopenmp susan_omp.c -lm -o susan_omp

pth: susan_pth.c
	gcc -o susan_pth susan_pth.c -lpthread

clean:
	rm susan_serial susan_mpi susan_omp
