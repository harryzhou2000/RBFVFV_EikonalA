

LIB=-lmkl_sequential -lmkl_core -lmkl_intel_lp64

# OPT=-g
OPT=-O3

OBJ=FieldSolver.o Grid.o Math.o Parameter.o


main.exe: main.cpp $(OBJ)
	icpc $^ -o $@ $(OPT) $(LIB)


%.o: %.cpp
	icpc $^ -c -o $@ $(OPT)

.PHONY: clean
clean:
	rm -rf *.o *.exe
