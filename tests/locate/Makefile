
INCS = -I../../include 
LIBS = 
CFLAGS = -g -O2
LDFLAGS = -L../../bin 
OBJS = locate.o
TARGET = ../../bin/locate

all: $(OBJS)
	g++ $(LDFLAGS) $(LIBS) $(OBJS) -o $(TARGET) 

.cpp.o:
	g++ $(INCS) $(CFLAGS) -c $< -o $@

locate.o: ../../include/compressed_quadtree.h ../../include/skip_quadtree.h

clean:
	rm *.o $(TARGET) 
