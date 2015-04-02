MPICPP=$(MPI_BIN)/mpic++
LDFLAGS=$(LD_FLAGS) -lfftw3_mpi -lfftw3 -lm
CFLAGS=-c -Wall

all: main

main: main.o
	$(MPICPP) $(LDFLAGS) main.o -o main

main.o: main.cpp
	$(MPICPP) $(CFLAGS) $(LDFLAGS) main.cpp 

tests: testing.o
	$(MPICPP) testing.o -o testing

testing.o: testing.cpp
	$(MPICPP) $(CFLAGS) testing.cpp

clean:
	rm *.o main test
