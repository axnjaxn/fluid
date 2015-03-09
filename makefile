CFLAGS = `byteimage-config --cflags` -O3 -Wno-unused-result -std=c++11

all: Simulate

fluidsim.o: fluidsim.h fluidsim.cpp
	$(CXX) fluidsim.cpp -c $(CFLAGS) 

main.o: fluidsim.h main.cpp
	$(CXX) main.cpp -c $(CFLAGS)

Simulate: fluidsim.o main.o
	$(CXX) fluidsim.o main.o -o Simulate `byteimage-config --libs` 

run: Simulate
	./Simulate

clean:
	rm -f *~ *.o Simulate
