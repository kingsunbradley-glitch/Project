CXX ?= g++
ROOT_CONFIG ?= root-config
CPPFLAGS += -Iinclude $(shell $(ROOT_CONFIG) --cflags)
CXXFLAGS += -O2 -Wall -Wextra -Wpedantic
LDLIBS += $(shell $(ROOT_CONFIG) --libs)
SOURCES := $(wildcard src/*.cpp)
wave_reader: $(SOURCES) $(wildcard include/*.h)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(SOURCES) $(LDLIBS) -o $@
.PHONY: all test clean
all: wave_reader
test: wave_reader
	python3 tests/validate.py
clean:
	rm -f wave_reader
