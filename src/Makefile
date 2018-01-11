CXXFLAGS = -O3 -fPIC 
COMPILER = gcc
objects = dwt.o 

libfft.so : $(objects)
	$(COMPILER) $(CXXFLAGS) -shared -o libdwt.so $(objects) -lm

dwt.o: dwt.c dwt.h
	$(COMPILER) $(CXXFLAGS) -c dwt.c 

clean: 
	rm *.o *.so