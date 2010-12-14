
INCS = -I
LIBS = -lann
CFLAGS = -g -O2
LDFLAGS = -g
OBJS = knn.o
TARGET = ../../bin/ann-knn

all: $(OBJS)
	g++ $(LDFLAGS) $(LIBS) $(OBJS) -o $(TARGET) 

.cpp.o:
	g++ $(INCS) $(CFLAGS) -c $< -o $@

clean:
	rm *.o $(TARGET) 
