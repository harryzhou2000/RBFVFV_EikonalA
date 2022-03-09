

LIB=-lmkl_sequential -lmkl_core -lmkl_intel_lp64

# OPT=-g
OPT= -O3 -g

OBJ=FieldSolver.o Grid.o Math.o Parameter.o

CPC=icpc

main.exe: main.cpp $(OBJ)
	$(CPC) $^ -o $@ $(OPT) $(LIB)

test1.exe: test1.cpp
	$(CPC) $^ -o $@

test2.exe: test2.cpp Math.o
	$(CPC) $^ -o $@ $(LIB)

%.o: %.cpp
	$(CPC) $^ -c -o $@ $(OPT)

.PHONY: clean
clean:
	rm -f *.o *.exe
