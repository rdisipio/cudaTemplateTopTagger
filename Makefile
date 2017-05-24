CC = gcc -std=c++11
DEBUG = -g
CUDA_DIR=$(SCINET_CUDA_BASE)
CUDA_INCDIR=$(SCINET_CUDA_INC)
CUDA_LIBDIR=$(SCINET_CUDA_LIB)

#INCFLAG=-I$(CUDA_DIR)/include -I.. -I./cub-1.4.1
CUDA_INCS=-I$(CUDA_INCDIR)

# CUDA code generation flags
#GENCODE_FLAGS   := -gencode arch=compute_35,code=sm_35
GENCODE_FLAGS = -arch=sm_35

NVCC = nvcc
NVCCFLAGS = $(DEBUG) -m64 $(INCFLAG) -Xcompiler '-fPIC' -std=c++11
#NVCCFLAGS = -g -m64 -Xcompiler '-fPIC'

CUDA_LIBS = -L$(CUDA_LIBDIR) -lcuda -lcudart -lcudadevrt -lcurand
#-lcublas -lcublas_device 
#CUDA_LIBS = -L$(CUDA_LIBDIR) -lcuda -lcudart -lcudadevrt 

#ROOT_FLAGS=$(shell root-config --cflags)
ROOT_INCS=-I$(shell root-config --incdir)
#ROOT_INCS=-I/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.02.12-x86_64-slc6-gcc48-opt/include
#ROOT_FLAGS=-std=c++11 -m64 -I/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.02.12-x86_64-slc6-gcc48-opt/include
ROOT_LIBS=-L/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.10-x86_64-slc6-gcc48-opt/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lm -ldl
#ROOT_LIBS=$(shell root-config --libs)

all: exec

%.o: %.C
	$(CC) $(DEBUG) -fPIC $(ROOT_INCS)  $(INCFLAG) -c $< -o $@

%.o: %.cxx
	$(CC) $(DEBUG) -fPIC $(ROOT_FLAGS) $(FASTJET_FLAGS) $(INCFLAG) $(CUDA_INCS) -c $< -o $@

%.o: %.cu
	$(NVCC) $(GENCODE_FLAGS) $(NVCCFLAGS) $(CUDA_INCS) $(ROOT_INCS) -lineinfo -dc $< -o $@

exec: run_gTTT.o
	$(NVCC) $(GENCODE_FLAGS) $(ROOT_LIBS) $(CUDA_LIBS) $< -o run_gTTT

clean:
	@rm -f *.o *~
