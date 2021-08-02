TARGETS = test_random

all: ${TARGETS}

clean:
	rm -f *.o ${TARGETS}

test_random:test_random.cpp SolveDALM_mex.cpp
	g++ test_random.cpp SolveDALM.h SolveDALM_mex.cpp -DBLAS_IMPLEMENTATION_MATLAB -I/Applications/MATLAB_R2009b.app/extern/include/ -L/Applications/MATLAB_R2009b.app/bin/maci64 -I../../util -lmwblas -lmwlapack -lm -o test_random