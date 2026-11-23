# =======================================================================
#   DSSD Calibration Pipeline Makefile
# =======================================================================

# --- 1. 编译器设置 ---
CXX       := g++
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --glibs) -lROOTDataFrame -lRooFit -lSpectrum

CXXFLAGS  := -O2 -Wall -fPIC -pthread $(ROOTCFLAGS)
LDFLAGS   := $(ROOTLIBS)

# --- 2. 目标文件名 ---
TARGET_PRE   := DSSD_Pre
TARGET_NORM  := DSSD_Norm
TARGET_CALIB := DSSD_Calib

SRC_PRE      := Preselect_Main.cpp
SRC_NORM     := Normalize_Main.cpp
SRC_CALIB    := Calibrator_Main.cpp Calibrator.cpp
OBJ_CALIB    := Calibrator.o

# --- 3. 编译规则 ---

all: $(TARGET_PRE) $(TARGET_NORM) $(TARGET_CALIB)
	@echo "-------------------------------------------"
	@echo "Build Success! Executables created:"
	@echo "  0. $(TARGET_PRE)   (Pre-selection)"
	@echo "  1. $(TARGET_NORM)  (Normalization)"
	@echo "  2. $(TARGET_CALIB) (Calibration)"
	@echo "-------------------------------------------"

$(TARGET_PRE): Preselect_Main.cpp Config.h
	@echo "[Compiling Pre] $@"
	$(CXX) $(CXXFLAGS) -o $@ Preselect_Main.cpp $(LDFLAGS)

$(TARGET_NORM): Normalize_Main.cpp Config.h
	@echo "[Compiling Norm] $@"
	$(CXX) $(CXXFLAGS) -o $@ Normalize_Main.cpp $(LDFLAGS)

$(OBJ_CALIB): Calibrator.cpp Calibrator.h Config.h
	@echo "[Compiling Object] $@"
	$(CXX) $(CXXFLAGS) -c Calibrator.cpp -o $@

$(TARGET_CALIB): Calibrator_Main.cpp $(OBJ_CALIB) Config.h
	@echo "[Compiling Cali] $@"
	$(CXX) $(CXXFLAGS) -o $@ Calibrator_Main.cpp $(OBJ_CALIB) $(LDFLAGS)


# --- 4. 运行指令 ---

.PHONY: Pre Norm Cali run clean help clean_Pre clean_Norm clean_Cali clean_all

Pre: $(TARGET_PRE)
	@echo "\n>>> Running Step 0: Pre-selection..."
	./$(TARGET_PRE)

Norm: $(TARGET_NORM)
	@echo "\n>>> Running Step 1: Normalization..."
	./$(TARGET_NORM)

Cali: $(TARGET_CALIB)
	@echo "\n>>> Running Step 2: Recalibration..."
	./$(TARGET_CALIB)

run: $(TARGET_PRE) $(TARGET_NORM) $(TARGET_CALIB)
	@echo "\n>>> Starting Full Pipeline..."
	@echo ">>> [1/3] Executing Pre"
	./$(TARGET_PRE)
	@echo "\n>>> [2/3] Executing Norm"
	./$(TARGET_NORM)
	@echo "\n>>> [3/3] Executing Cali"
	./$(TARGET_CALIB)
	@echo "\n>>> Pipeline Finished Successfully."

# --- 5. 清理规则 (Detailed Clean) ---

# 清除 Preselect 产生的中间文件
clean_Pre:
	@echo "Cleaning Pre-selection outputs..."
	rm -f SS032_decay_only.root SS032_implant_only.root

# 清除 Normalize 产生的中间文件
clean_Norm:
	@echo "Cleaning Normalization outputs..."
	rm -f SS032_Normalize_Params.txt Diagnose_Normalize.root

# 清除 Calibrator 产生的中间文件
clean_Cali:
	@echo "Cleaning Calibration outputs..."
	rm -f ener_Recal.dat Diagnose_Calibration.root peaks_list.txt

# 清除所有可执行文件和 .o 文件
clean:
	@echo "Cleaning executables and objects..."
	rm -f $(TARGET_PRE) $(TARGET_NORM) $(TARGET_CALIB) *.o

# 清除一切（包括所有生成的数据文件）
clean_all: clean clean_Pre clean_Norm clean_Cali
	@echo "All generated files have been removed."

help:
	@echo "Available commands:"
	@echo "  make         - Compile all programs"
	@echo "  make Pre     - Run Step 0 (Pre-selection)"
	@echo "  make Norm    - Run Step 1 (Normalization)"
	@echo "  make Cali    - Run Step 2 (Calibration)"
	@echo "  make run     - Run FULL sequence"
	@echo "  make clean   - Remove executables only"
	@echo "  make clean_Pre  - Remove Pre-select data files"
	@echo "  make clean_Norm - Remove Normalize data files"
	@echo "  make clean_Cali - Remove Calibrate data files"
	@echo "  make clean_all  - Remove EVERYTHING (Reset)"