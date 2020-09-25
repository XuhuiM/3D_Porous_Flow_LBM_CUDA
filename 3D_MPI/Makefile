CC=mpicc
CFLAGS=-c -O3
SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=flow

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $^ -o $@
.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm -rf *.o  *.dat flow


