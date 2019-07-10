# include folders
BOOST_INCLUDE_DIR = 
# lib folders
BOOST_LIB_DIR = 

CC = gcc
CXX = g++
CXXFLAGS = -std=c++11 $(INCLUDES) -g
INCLUDES = -I $(BOOST_INCLUDE_DIR) 
LDLIBS = -lz 

SRCS_CATE = src/CategorizeSimulation.cpp src/Transcript.cpp

all: bin/categorizesimulation

bin/categorizesimulation: $(subst .cpp,.o,$(SRCS_CATE))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDLIBS) 

clean:
	rm -f bin/* src/*.o
