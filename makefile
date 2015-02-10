all: app

app: main.cpp
	$(CXX) main.cpp -o app `byteimage-config --cflags --libs` -O3 -Wno-unused-result

run: app
	./app

clean:
	rm -f *~ *.o app
