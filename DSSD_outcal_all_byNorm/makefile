# ==========================
#  SS032 3-peak calibration
#  build:
#    make Pres | Norm | Cali | all
#  run (build + run):
#    make runPres  MAP=run00001_map.root PRE=preselected.root
#    make runNorm  PRE=preselected.root NORM=SS032_Normalize_Params.txt NDIAG=Diagnose_Normalize.root
#    make runCali  PRE=preselected.root NORM=SS032_Normalize_Params.txt OUT=ener_cal.dat DMAIN=diagMain.root DED=diaged.root
# ==========================

CXX      := g++
CXXFLAGS := -O2 -std=c++17 -Wall -Wextra

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)
SPECTRUM_LIB := -lSpectrum

PRES_BIN := Preselect_Main
NORM_BIN := Normalize_Main
CALI_BIN := Calibrator_Main

PRES_SRC := Preselect_Main.cpp
NORM_SRC := Normalize_Main.cpp
CALI_SRC := Calibrator_Main.cpp Calibrator.cpp

.PHONY: all Pres Norm Cali clean runPres runNorm runCali

all: Pres Norm Cali

Pres: $(PRES_BIN)
Norm: $(NORM_BIN)
Cali: $(CALI_BIN)

$(PRES_BIN): $(PRES_SRC)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $^ $(ROOTLIBS) -o $@

$(NORM_BIN): $(NORM_SRC)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $^ $(ROOTLIBS) -o $@

$(CALI_BIN): $(CALI_SRC)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $^ $(ROOTLIBS) $(SPECTRUM_LIB) -o $@

clean:
	rm -f $(PRES_BIN) $(NORM_BIN) $(CALI_BIN) *.o *.so *.d

# --------------------------
# Run targets (build + run)
# --------------------------

# Defaults (override on command line)
MAP   ?= run00001_map.root
PRE   ?= preselected.root

NORM  ?= SS032_Normalize_Params.txt
NDIAG ?= Diagnose_Normalize.root

OUT   ?= ener_cal.dat
DMAIN ?= diagMain.root
DED   ?= diaged.root

runPres: Pres
	./$(PRES_BIN) "$(MAP)" "$(PRE)"

runNorm: Norm
	./$(NORM_BIN) "$(PRE)" "$(NORM)" "$(NDIAG)"

runCali: Cali
	./$(CALI_BIN) "$(PRE)" "$(NORM)" "$(OUT)" "$(DMAIN)" "$(DED)"
