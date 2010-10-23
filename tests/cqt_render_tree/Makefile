
INCS = -I../../include 
LIBS = 
CFLAGS = 
LDFLAGS = 
OBJS = render_tree.o
TARGET = ../../bin/render-tree

all: $(OBJS)
	g++ $(LDFLAGS) $(LIBS) $(OBJS) -o $(TARGET) 

.cpp.o:
	g++ $(INCS) $(CFLAGS) -c $< -o $@

clean:
	rm *.o $(TARGET) 
