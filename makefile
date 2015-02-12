all: Simulate

fluidsim.o: fluidsim.h fluidsim.cpp
	$(CXX) fluidsim.cpp -c `byteimage-config --cflags` -O3 -Wno-unused-result

main.o: fluidsim.h main.cpp
	$(CXX) main.cpp -c `byteimage-config --cflags` -O3 -Wno-unused-result

Simulate: fluidsim.o main.o
	$(CXX) fluidsim.o main.o -o Simulate `byteimage-config --libs` 

run: Simulate
	./Simulate

clean:
	rm -f *~ *.o Simulate
