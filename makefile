all: floydWarshall

clean:
	rm *.o a.out
	
floydWarshall: floydWarshall.o
	mpicc floydWarshall.o -o floydWarshall

floydWarshall.o: floydWarshall.c
	mpicc -c floydWarshall.c