# include folders
BOOST_INCLUDE_DIR = 
# lib folders
BOOST_LIB_DIR = 

CC = gcc
CXX = g++
CXXFLAGS = -std=c++11 $(INCLUDES) -g
INCLUDES = -I $(BOOST_INCLUDE_DIR) 
LDLIBS = -lz 

SRCS_ASEM = src/AssemblyPost.cpp
SRCS_CATE = src/CategorizeSimulation.cpp src/Transcript.cpp

all: bin/assemblypost bin/categorizesimulation

bin/assemblypost: $(subst .cpp,.o,$(SRCS_ASEM))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS) -Wl,-rpath,$(RPATH)

bin/categorizesimulation: $(subst .cpp,.o,$(SRCS_CATE))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDLIBS) 

clean:
	rm -f bin/* src/*.o
