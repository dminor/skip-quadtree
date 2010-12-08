
INCS = -I../../include 
LIBS = 
CFLAGS = -g
LDFLAGS = -L../../bin 
OBJS = knn.o
TARGET = ../../bin/sqt-knn

all: $(OBJS)
	g++ $(LDFLAGS) $(LIBS) $(OBJS) -o $(TARGET) 

.cpp.o:
	g++ $(INCS) $(CFLAGS) -c $< -o $@

knn.o: ../../include/skip_quadtree.h

clean:
	rm *.o $(TARGET) 
