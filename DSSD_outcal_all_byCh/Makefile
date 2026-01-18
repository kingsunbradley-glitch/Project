CXX := g++
CXXFLAGS := -O2 -std=c++17 -Wall -Wextra
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

BIN := PerChannelCalibrator
SRC := PerChannelCalibrator.cpp

.PHONY: all clean run

all: $(BIN)

$(BIN): $(SRC)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $^ $(ROOTLIBS) -lSpectrum -o $@

clean:
	rm -f $(BIN) *.o

# Defaults (override on command line)
IN   ?= preselected.root
TREE ?= tr_map
ORIG ?= originalSpectrum.root
CAL  ?= calibratedSpectrum.root
DAT  ?= ener_cal.dat

run: $(BIN)
	./$(BIN) "$(IN)" "$(TREE)" "$(ORIG)" "$(CAL)" "$(DAT)"
