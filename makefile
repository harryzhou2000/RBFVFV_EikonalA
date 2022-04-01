all: main.exe



#### Intel Compiler: 
### doesn't work now
# LIB=-lmkl_sequential -lmkl_core -lmkl_intel_lp64 -fopenmp
# FPFLAGS=-fp-model precise -inline-forceinline
# CPC=icpc

#### GCC or CLANG/LLVM: 
FPFLAGS=-DUSE_OPEN_BLAS -fopenmp 
CPC=/usr/bin/clang++ 
# CPC=/usr/bin/g++
LIB=-lopenblas -llapacke



OPT=-g
OPT= -O3 -g 
OPT= -O3 -DNDEBUG # NDEBUG macro for c assert library

OBJ:=FieldSolver.o Grid.o Math.o Parameter.o main.o

# DEP:=$(patsubst %, %.d, $(OBJ))
DEP:=$(OBJ:.o=.d) # main.d, Math.d etc

OBJ:=$(OBJ) MathHard.o ## Header dependencies not tracked for MathHard.o, take care if it's affected

-include $(DEP) ## '-' is for first compliling that has not (enough) .d files


main.exe:  $(OBJ)
	$(CPC) $^ -o $@ $(OPT) $(LIB) $(FPFLAGS)

test1.exe: test1.cpp
	$(CPC) $^ -o $@

test2.exe: test2.cpp Math.o MathHard.o
	$(CPC) $^ -o $@ $(LIB) $(FPFLAGS)


%.o: %.cpp 
# mind that only first input is compiled for other dependencies are included files
# mind that -MMD instead of -MM to actually compile it
	$(CPC) $< -c -o $@ $(OPT)  $(FPFLAGS) -MMD 

# %.d: %.cpp
# 	%(CPC) -MM $^ > $@ $(OPT)  $(FPFLAGS)

# make A to see the order of dependencies
C:

B:

D:

M:

A:  D M 

A : C

A : B
	echo $(^) 


.PHONY: clean
clean:
	rm -f *.o *.d *.exe
