



#### Intel Compiler:
# LIB=-lmkl_sequential -lmkl_core -lmkl_intel_lp64 -fopenmp
# FPFLAGS=-fp-model precise -inline-forceinline
# CPC=icpc

#### GCC: 
FPFLAGS=
CPC=g++
LIB=-fopenmp -DUSE_OPEN_BLAS -lopenblas -llapacke



OPT=-g
OPT= -O3 -g 
# OPT= -O3 -DNDEBUG

OBJ=FieldSolver.o Grid.o Math.o Parameter.o



main.exe: main.cpp $(OBJ)
	$(CPC) $^ -o $@ $(OPT) $(LIB) $(FPFLAGS)

test1.exe: test1.cpp
	$(CPC) $^ -o $@

test2.exe: test2.cpp Math.o
	$(CPC) $^ -o $@ $(LIB) $(FPFLAGS)

%.o: %.cpp
	$(CPC) $^ -c -o $@ $(OPT) $(LIB) $(FPFLAGS)

.PHONY: clean
clean:
	rm -f *.o *.exe
