# ===================================================================
# Makefile for DSSD Data Processing Workflow
# ===================================================================

# 编译器和选项
CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17

# ROOT 库和编译选项
ROOT_CFLAGS = $(shell root-config --cflags)
ROOT_LIBS = $(shell root-config --libs)

# RDataFrame 特定的库
RDF_LIB = -lROOTDataFrame -lTreePlayer
# 最终的编译和链接选项
FINAL_CXXFLAGS = $(CXXFLAGS) $(ROOT_CFLAGS)
FINAL_LIBS = $(ROOT_LIBS) $(RDF_LIB)

# 定义要生成的可执行文件
EXECUTABLES = DSSD_normalize sorter

# --- 从 config.h 读取清理所需的文件和目录 ---
CONFIG_H = config.h
PARAM_FILE = $(shell grep 'normalize_param_file' $(CONFIG_H) | sed -e 's/.*"\(.*\)"[^"]*$$/\1/')
PLOT_DIR = $(shell grep 'output_plot_dir' $(CONFIG_H) | sed -e 's/.*"\(.*\)"[^"]*$$/\1/')

# --- 规则 ---

# "all" 是默认目标, 编译所有可执行文件
all: $(EXECUTABLES)

# 编译 DSSD_normalize 的规则
DSSD_normalize: DSSD_normalize.cpp $(CONFIG_H)
	$(CXX) $(FINAL_CXXFLAGS) $< -o $@ $(FINAL_LIBS)

# 编译 sorter 的规则
sorter: sorter.cpp $(CONFIG_H)
	$(CXX) $(FINAL_CXXFLAGS) $< -o $@ $(FINAL_LIBS)

# "run" 规则：一键执行整个分析流程
run: all
	@echo "--- [STEP 1/2] Running Pre-sorterer to create snapshot file... ---"
	./sorter
	@echo "\n--- [STEP 2/2] Running Normalizer on the snapshot... ---"
	./DSSD_normalize
	@echo "\n--- DSSD Normalization Workflow Finished! ---"

# "clean" 规则：删除所有生成的文件
.PHONY: clean all run
clean:
	@echo "Cleaning up compiled executables..."
	rm -f $(EXECUTABLES) *.o
	@echo "Cleaning up intermediate and output files..."
	rm -f sorted.root
	rm -f $(PARAM_FILE)
	rm -rf $(PLOT_DIR)
	@echo "Cleanup complete."

