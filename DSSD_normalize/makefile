# ===================================================================
# Makefile for DSSD Normalization Program
# ===================================================================

# 编译器和选项
CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17

# ROOT 库和编译选项
ROOT_CFLAGS = $(shell root-config --cflags)
ROOT_LIBS = $(shell root-config --libs)

# RDataFrame 特定的库
RDF_LIB =  -lTreePlayer
#-lRDataFrame
# 最终的编译和链接选项
FINAL_CXXFLAGS = $(CXXFLAGS) $(ROOT_CFLAGS)
FINAL_LIBS = $(ROOT_LIBS) $(RDF_LIB)

# 定义要生成的可执行文件
EXECUTABLES = DSSD_normalize

# --- 自动从 config_dssd.h 读取输出文件和目录路径 ---
CONFIG_H = config.h
PARAM_FILE = $(shell grep 'normalize_param_file' $(CONFIG_H) | sed -e 's/.*"\(.*\)"[^"]*$$/\1/')
PLOT_DIR = $(shell grep 'output_plot_dir' $(CONFIG_H) | sed -e 's/.*"\(.*\)"[^"]*$$/\1/')
EXPERIMENT_ID = $(shell grep 'EXPERIMENT_ID' $(CONFIG_H) | sed -e 's/.*"\(.*\)"[^"]*$$/\1/')
DIAGNOSTIC_ROOT_FILE = $(PLOT_DIR)/$(EXPERIMENT_ID)_Normalize_Diagnostic.root


# --- 规则 ---

# "all" 是默认目标
all: $(EXECUTABLES)

# 编译 DSSD_normalize 的规则
DSSD_normalize: DSSD_normalize.cpp $(CONFIG_H)
	$(CXX) $(FINAL_CXXFLAGS) $< -o $@ $(FINAL_LIBS)

# "run" 规则：先编译，然后执行
run: all
	@echo "--- Running DSSD Normalization: Generating k/b parameters and diagnostics... ---"
	./DSSD_normalize
	@echo "\n--- DSSD Normalization Workflow Finished! ---"

# "clean" 规则：删除所有生成的文件
.PHONY: clean all run
clean:
	@echo "Cleaning up compiled executables and object files..."
	rm -f $(EXECUTABLES) *.o
	@echo "Cleaning up generated data files and plots defined in $(CONFIG_H)..."
	rm -f $(PARAM_FILE)
	rm -f $(DIAGNOSTIC_ROOT_FILE)
	rm -rf $(PLOT_DIR)
	@echo "Cleanup complete."