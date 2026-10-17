CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -pedantic -fopenmp
LDLIBS ?= -lstdc++fs -fopenmp

.PHONY: all test clean gpu

all: maree2006

maree2006: maree2006.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDLIBS)

test: maree2006
	./maree2006 --self-test
	python3 test_model.py

clean:
	rm -f maree2006 maree2006_gpu

gpu: maree2006_gpu

maree2006_gpu: maree2006_gpu.cu
	nvcc -O3 -std=c++17 -arch=sm_80 --use_fast_math -lineinfo $< -o $@ -lstdc++fs
