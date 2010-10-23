
INCS = -I../../include 
LIBS = 
CFLAGS = -g
LDFLAGS = -L../../bin 
OBJS = knn.o
TARGET = ../../bin/cqt_knn

all: $(OBJS)
	g++ $(LDFLAGS) $(LIBS) $(OBJS) -o $(TARGET) 

.cpp.o:
	g++ $(INCS) $(CFLAGS) -c $< -o $@

clean:
	rm *.o $(TARGET) 
