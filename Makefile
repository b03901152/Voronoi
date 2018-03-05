CC = g++-7
#DEBUG = -O3 -DNDEBUG
DEBUG = -O3 
CFLAGS = -lGL -lGLU -lglut -g -std=c++11 -Wall $(DEBUG)

all:
	$(CC) src/main.cpp  $(CFLAGS)  -o Voronoi
clean:
	rm -rf *.o *.gprof Voronoi lib/*.a lib/*.o bin/*