all: app

fluidsim.o: fluidsim.h fluidsim.cpp
	$(CXX) fluidsim.cpp -c `byteimage-config --cflags` -O3 -Wno-unused-result

main.o: fluidsim.h main.cpp
	$(CXX) main.cpp -c `byteimage-config --cflags` -O3 -Wno-unused-result

app: fluidsim.o main.o
	$(CXX) fluidsim.o main.o -o app `byteimage-config --libs` 

run: app
	./app

clean:
	rm -f *~ *.o app
