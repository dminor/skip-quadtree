
INCS = -I../../../kdtree/include 
LIBS = 
CFLAGS = -g -O2
LDFLAGS = -L../../bin 
OBJS = knn.o
TARGET = ../../bin/kdtree-knn

all: $(OBJS)
	g++ $(LDFLAGS) $(LIBS) $(OBJS) -o $(TARGET) 

.cpp.o:
	g++ $(INCS) $(CFLAGS) -c $< -o $@

knn.o: ../../../kdtree/include/kdtree.h

clean:
	rm *.o $(TARGET) 
