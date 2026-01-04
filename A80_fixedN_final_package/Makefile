# 简易 Makefile（ROOT 工程常用写法）
CXX := g++
CXXFLAGS := -O2 -std=c++17 $(shell root-config --cflags)
LDLIBS := $(shell root-config --libs)

SRCS := main.cpp RunProcessor.cpp PeakFinder.cpp Metrics.cpp QaIO.cpp
OBJS := $(SRCS:.cpp=.o)

all: process_runs_A80

process_runs_A80: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) process_runs_A80

.PHONY: all clean
