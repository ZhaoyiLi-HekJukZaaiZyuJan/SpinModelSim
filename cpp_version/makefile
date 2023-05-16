DIRS := . ../src
SOURCES := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cpp))
OBJS := $(patsubst %.cpp, %.o, $(SOURCES))

CFLAGS := -std=c++17 -no-pie -fmax-errors=5 -O3\
	-I ../src \
	-I ../src/libs/
CXX ?= g++
INCLUDES :=
LIBDIR :=

all: simulate

simulate: ${OBJS}
	$(CXX) $(CFLAGS) -o $@ ${OBJS} ${LIBS} \
	

.cpp.o:
	$(CXX) $(CFLAGS) ${INCLUDES} $< -c -o $@ | echo $<

clean:
	rm -f ${OBJS} simulate