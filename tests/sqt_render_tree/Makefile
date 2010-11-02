
INCS = -I../../include 
LIBS = 
CFLAGS = -g
LDFLAGS = 
OBJS = render_tree.o
TARGET = ../../bin/sqt-render-tree

all: $(OBJS)
	g++ $(LDFLAGS) $(LIBS) $(OBJS) -o $(TARGET) 

.cpp.o:
	g++ $(INCS) $(CFLAGS) -c $< -o $@

render_tree.o: ../../include/skip_quadtree.h

clean:
	rm *.o $(TARGET) 
