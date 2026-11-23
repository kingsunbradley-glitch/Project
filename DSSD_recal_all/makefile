# =======================================================================
#   DSSD Calibration Pipeline Makefile
# =======================================================================

# 显式指定默认目标。必须写：本文件第一个规则是 .setup_state（以点开头，
# GNU make 不把它当默认目标），再往下的规则是 FORCE 的占位规则，不加这一行
# 时裸敲 make 会去打 FORCE，输出 "Nothing to be done for 'FORCE'" 且什么也不编译。
.DEFAULT_GOAL := all

# --- 1. 编译器设置 ---
CXX       := g++
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --glibs) -lROOTDataFrame -lSpectrum

CXXFLAGS  := -O2 -Wall -fPIC -pthread $(ROOTCFLAGS)
LDFLAGS   := $(ROOTLIBS)

# --- 1b. 与 ANA_SHANS2 共享的 setup.h ---
# 本程序的道号区间(DSSD_X/Y/YH)和最终刻度文件名不再自己写死，而是从 ANA 那套
# 共享 setup.h 推导(见 Config.h 第 0 节)。这里只负责把它所在的目录塞进搜索路径；
# 目录不存在或没有 setup.h 时，Config.h 会退回内置兜底值并在启动时打印 [INFO]。
# 换机器/换版本改这一行即可：make ANA_SETUP_DIR=/path/to/Cal_Auto_DSSD_Normalize
ANA_SETUP_DIR ?= /home/evalie/Project/Function/ANA/ANA_SHANS2_202602_v1.07.4/Cal_Auto_All_202412/Cal_Auto_DSSD_Normalize
CXXFLAGS  += -I$(ANA_SETUP_DIR)

# --- 1c. setup.h 的两种用法 ---
#  方式 A (默认)：只靠上面的 -I 直接读共享文件 —— 全机只有一份定义，上游一改
#                立刻生效，不会出现"本地副本过期"的坑。
#  方式 B：make sync-setup —— 把共享 setup.h 快照到本目录，与 ANA 其它程序
#                "每个程序目录各带一份"的约定一致。带引号的 #include "setup.h"
#               会先找本目录，所以本地副本优先于 -I。
#  撤销方式 B：make unsync-setup
SETUP_LOCAL := setup.h
SETUP_STAMP := setup_stamp.txt
SETUP_STATE := .setup_state

# 实际生效的是哪一份：有本地快照就用它，否则用 -I 指向的那份
SETUP_MODE  := $(if $(wildcard $(SETUP_LOCAL)),local,include)
ifeq ($(SETUP_MODE),local)
  SETUP_EFFECTIVE := $(SETUP_LOCAL)
else
  SETUP_EFFECTIVE := $(wildcard $(ANA_SETUP_DIR)/setup.h)
endif

# 状态戳：内容 = 模式 + 实际生效 setup.h 的 md5。内容变了才更新 mtime，
# 从而让下面的 .o/可执行文件重编。这样"换模式"和"上游 setup.h 改了"都不会
# 静默沿用旧目标文件。
SETUP_SIG   := $(SETUP_MODE):$(shell [ -f "$(SETUP_EFFECTIVE)" ] && md5sum "$(SETUP_EFFECTIVE)" | cut -d' ' -f1 || echo none)

$(SETUP_STATE): FORCE
	@echo '$(SETUP_SIG)' > $@.tmp
	@cmp -s $@.tmp $@ || { mv $@.tmp $@; echo "[setup.h] 配置变化 -> $(SETUP_SIG)"; }
	@rm -f $@.tmp
FORCE:

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
	@if [ -f $(SETUP_LOCAL) ]; then \
	    echo "Setup.h : 本地快照 ./$(SETUP_LOCAL) (make unsync-setup 可撤销)"; \
	elif [ -f "$(ANA_SETUP_DIR)/setup.h" ]; then \
	    echo "Setup.h : -I $(ANA_SETUP_DIR)"; \
	else \
	    echo "Setup.h : [WARN] 找不到 $(ANA_SETUP_DIR)/setup.h，使用 Config.h 内置兜底道号"; fi
	@echo "-------------------------------------------"
	@echo "Build Success! Executables created:"
	@echo "  0. $(TARGET_PRE)   (Pre-selection)"
	@echo "  1. $(TARGET_NORM)  (Normalization)"
	@echo "  2. $(TARGET_CALIB) (Calibration)"
	@echo "-------------------------------------------"

$(TARGET_PRE): Preselect_Main.cpp Config.h $(SETUP_STATE)
	@echo "[Compiling Pre] $@"
	$(CXX) $(CXXFLAGS) -o $@ Preselect_Main.cpp $(LDFLAGS)

$(TARGET_NORM): Normalize_Main.cpp Config.h $(SETUP_STATE)
	@echo "[Compiling Norm] $@"
	$(CXX) $(CXXFLAGS) -o $@ Normalize_Main.cpp $(LDFLAGS)

$(OBJ_CALIB): Calibrator.cpp Calibrator.h Config.h $(SETUP_STATE)
	@echo "[Compiling Object] $@"
	$(CXX) $(CXXFLAGS) -c Calibrator.cpp -o $@

$(TARGET_CALIB): Calibrator_Main.cpp $(OBJ_CALIB) Config.h $(SETUP_STATE)
	@echo "[Compiling Cali] $@"
	$(CXX) $(CXXFLAGS) -o $@ Calibrator_Main.cpp $(OBJ_CALIB) $(LDFLAGS)


# --- 4. 运行指令 ---

.PHONY: FORCE Pre Norm Cali run clean help clean_Pre clean_Norm clean_Cali clean_all \
        sync-setup unsync-setup

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

# --- 4b. setup.h 快照管理 ---

sync-setup:
	@if [ ! -f "$(ANA_SETUP_DIR)/setup.h" ]; then \
	    echo "[ERROR] 找不到 $(ANA_SETUP_DIR)/setup.h"; \
	    echo "        用 make sync-setup ANA_SETUP_DIR=/path/to/Cal_Auto_DSSD_Normalize 指定源目录"; \
	    exit 1; fi
	@cp "$(ANA_SETUP_DIR)/setup.h" $(SETUP_LOCAL)
	@printf 'source : %s\nmd5    : %s\ncopied : %s\n' \
	    "$(ANA_SETUP_DIR)/setup.h" \
	    "$$(md5sum $(SETUP_LOCAL) | cut -d' ' -f1)" \
	    "$$(date '+%F %T %Z')" > $(SETUP_STAMP)
	@echo "--> 已把共享 setup.h 快照到 ./$(SETUP_LOCAL)："
	@sed 's/^/    /' $(SETUP_STAMP)
	@echo "--> 本地副本优先于 -I；上游 setup.h 更新后请重新 make sync-setup"
	@echo "--> 撤销：make unsync-setup"
	@rm -f $(SETUP_STATE)

unsync-setup:
	@rm -f $(SETUP_LOCAL) $(SETUP_STAMP)
	@echo "--> 已删除本地 setup.h 快照，恢复 -I 直接读取："
	@echo "    $(ANA_SETUP_DIR)"
	@rm -f $(SETUP_STATE)

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
	rm -f ener_cal_Recal.dat Diagnose_Calibration.root peaks_list.txt

# 清除所有可执行文件和 .o 文件
clean:
	@echo "Cleaning executables and objects..."
	rm -f $(TARGET_PRE) $(TARGET_NORM) $(TARGET_CALIB) *.o $(SETUP_STATE)

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
	@echo "  make sync-setup   - 把共享 setup.h 快照到本目录(与 ANA 其它程序一致)"
	@echo "  make unsync-setup - 删除本地快照，恢复 -I 直接读取共享文件"