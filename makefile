



#### Intel Compiler: 
### doesn't work now
# LIB=-lmkl_sequential -lmkl_core -lmkl_intel_lp64 -fopenmp
# FPFLAGS=-fp-model precise -inline-forceinline
# CPC=icpc

#### GCC or CLANG/LLVM: 
FPFLAGS=
CPC=/usr/bin/clang++
# CPC=/usr/bin/g++
LIB=-fopenmp -DUSE_OPEN_BLAS -lopenblas -llapacke



OPT=-g
OPT= -O3 -g 
# OPT= -O3 -DNDEBUG

OBJ:=FieldSolver.o Grid.o Math.o Parameter.o main.o

# DEP:=$(patsubst %, %.d, $(OBJ))
DEP:=$(OBJ:.o=.d)

OBJ:=$(OBJ) MathHard.o ## Header dependencies not tracked for MathHard.o, take care if it's affected

-include $(DEP)

all:
	echo $(DEP)

main.exe:  $(OBJ)
	$(CPC) $^ -o $@ $(OPT) $(LIB) $(FPFLAGS)

test1.exe: test1.cpp
	$(CPC) $^ -o $@

test2.exe: test2.cpp Math.o
	$(CPC) $^ -o $@ $(LIB) $(FPFLAGS)

%.o: %.cpp
	$(CPC) $< -c -o $@ $(OPT)  $(FPFLAGS) -MMD

# %.d: %.cpp
# 	%(CPC) -MM $^ > $@ $(OPT)  $(FPFLAGS)



.PHONY: clean
clean:
	rm -f *.o *.d *.exe
