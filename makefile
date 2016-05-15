all: SamplingModel calculateOverpasses
SamplingModel:
	gcc -std=c++0x -Wall  -pedantic -o samplingModel sampling_model1.cpp -lgsl -lgslcblas  -I/opt/anaconda/include/ /opt/anaconda/lib/libgdal.so

calculateOverpasses: 
	gcc -std=c++0x -Wall  -pedantic -o calculateOverpasses calculateOverpasses.cpp -I/opt/anaconda/include/ /opt/anaconda/lib/libgdal.so
