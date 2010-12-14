
INCS = -I../../include 
LIBS = 
CFLAGS = -g -O2
LDFLAGS = -L../../bin 
OBJS = knn.o
TARGET = ../../bin/cqt-knn

all: $(OBJS)
	g++ $(LDFLAGS) $(LIBS) $(OBJS) -o $(TARGET) 

.cpp.o:
	g++ $(INCS) $(CFLAGS) -c $< -o $@

knn.o: ../../include/compressed_quadtree.h

clean:
	rm *.o $(TARGET) 
