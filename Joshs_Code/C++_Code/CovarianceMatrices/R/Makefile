LD=g++
LDFLAGS=-L/usr/lib64 -lfftw3 -g -L/usr/lib64 -lgsl -lgslcblas -L/home/jsdillon/Libraries/FFTW/lib -L/home/jsdillon/Libraries/GSL/lib 
CXXFLAGS= -g
IFLAGS = -I/home/jsdillon/Libraries/GSL/include/ -I/home/jsdillon/Libraries/FFTW/include/

R : R.o ../../CommonClasses/Specs.o
	$(LD) -o $@ $^ $(LDFLAGS)

R.o : R.cpp ../../CommonClasses/Specs.h
	$(LD) $(IFLAGS) -c $^ $(LDFLAGS)

Specs.o : ../../CommonClasses/Specs.cpp ../../CommonClasses/Specs.h
	$(LD) $(IFLAGS) -c $^ $(LDFLAGS)
	
.dummy clean :
	rm -f R *.o
