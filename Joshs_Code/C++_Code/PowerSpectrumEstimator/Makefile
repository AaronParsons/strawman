LD=g++
LDFLAGS=-L/usr/lib64 -lfftw3 -g -L/usr/lib64 -lgsl -lgslcblas -L/nfs/blank/h4213/jsdillon/Libraries/FFTW/lib -L/nfs/blank/h4213/jsdillon/Libraries/GSL/lib 
CXXFLAGS= -g
IFLAGS = -I/nfs/blank/h4213/jsdillon/Libraries/GSL/include/ -I/nfs/blank/h4213/jsdillon/Libraries/FFTW/include/

powerSpectrumEstimator : powerSpectrumEstimator.o ../CommonClasses/CVector.o ../CommonClasses/RDMatrix.o ../CommonClasses/Toeplitz.o ../CommonClasses/Specs.o
	$(LD) $(IFLAGS) -o  $@ $^ $(LDFLAGS)

powerSpectrumEstimator.o : powerSpectrumEstimator.cpp ../CommonClasses/CVector.h ../CommonClasses/RDMatrix.h ../CommonClasses/Toeplitz.h ../CommonClasses/Specs.h
	$(LD) -c $(IFLAGS) $^ $(LDFLAGS)

../CommonClasses/RDMatrix.o : ../CommonClasses/RDMatrix.cpp ../CommonClasses/RDMatrix.h ../CommonClasses/CVector.h ../CommonClasses/Specs.h

../CommonClasses/Toeplitz.o : ../CommonClasses/Toeplitz.cpp ../CommonClasses/Toeplitz.h ../CommonClasses/CVector.h ../CommonClasses/Specs.h
	$(LD) -c $(IFLAGS) $^ $(LDFLAGS)

../CommonClasses/CVector.o : ../CommonClasses/CVector.cpp ../CommonClasses/CVector.h ../CommonClasses/Specs.h
	$(LD) -c $(IFLAGS) $^ $(LDFLAGS)

../CommonClasses/Specs.o : ../CommonClasses/Specs.cpp ../CommonClasses/Specs.h

.dummy clean :
	rm -f powerSpectrumEstimator *.o