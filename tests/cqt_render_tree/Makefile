
INCS = -I../../include 
LIBS = 
CFLAGS = 
LDFLAGS = 
OBJS = render_tree.o
TARGET = ../../bin/cqt-render-tree

all: $(OBJS)
	g++ $(LDFLAGS) $(LIBS) $(OBJS) -o $(TARGET) 

.cpp.o:
	g++ $(INCS) $(CFLAGS) -c $< -o $@

render_tree.o: ../../include/compressed_quadtree.h

clean:
	rm *.o $(TARGET) 
