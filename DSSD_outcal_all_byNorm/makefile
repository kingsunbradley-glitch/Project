CXX ?= g++

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

CXXFLAGS += -O2 -std=c++17 $(ROOTCFLAGS)
LDLIBS   += $(ROOTLIBS) -lSpectrum

# -------------------------
# Runtime variables (override at make time)
# e.g. make Pres IN=run00001_map.root
#      make all IN=run00001_map.root
# -------------------------
IN       ?= run00001_map.root
PRESEL   ?= preselected.root
NORMTXT  ?= SS032_Normalize_Params.txt
NORMDIAG ?= Diagnose_Normalize.root
OUTDAT   ?= ener_cal.dat
DIAGMAIN ?= diagMain.root
DIAGED   ?= diaged.root

.PHONY: all clean Pres Norm Cali build

all: Pres Norm Cali

# -------- Build only --------
Preselect_Main: Preselect_Main.cpp Config.h
	$(CXX) $(CXXFLAGS) $< $(LDLIBS) -o $@

Normalize_Main: Normalize_Main.cpp Config.h
	$(CXX) $(CXXFLAGS) $< $(LDLIBS) -o $@

Calibrator_Main: Calibrator_Main.cpp Calibrator.cpp Calibrator.h Config.h
	$(CXX) $(CXXFLAGS) Calibrator.cpp Calibrator_Main.cpp $(LDLIBS) -o $@

# -------- Build + Run (requested) --------
Pres: Preselect_Main
	./Preselect_Main "$(IN)" "$(PRESEL)"

Norm: Normalize_Main
	./Normalize_Main "$(PRESEL)" "$(NORMTXT)" "$(NORMDIAG)"

Cali: Calibrator_Main
	./Calibrator_Main "$(PRESEL)" "$(NORMTXT)" "$(OUTDAT)" "$(DIAGMAIN)" "$(DIAGED)"

clean:
	rm -f Preselect_Main Normalize_Main Calibrator_Main
	rm -f "$(PRESEL)" "$(NORMTXT)" "$(NORMDIAG)" "$(OUTDAT)" "$(DIAGMAIN)" "$(DIAGED)"
