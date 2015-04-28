all: floydWarshall

clean:
	rm *.o floydWarshall
	
floydWarshall: floydWarshall.o
	mpicc floydWarshall.o -o floydWarshall

floydWarshall.o: floydWarshall.c
	mpicc -c floydWarshall.c