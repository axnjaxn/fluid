all: app

app: main.cpp
	$(CXX) main.cpp -o app `byteimage-config --cflags --libs`

run: app
	./app

clean:
	rm -f *~ *.o app
