# ===================================================================
# Makefile for SSD Calibration Programs (V2.0 - Enhanced Clean)
# ===================================================================

# 编译器和选项
CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++11

# ROOT 库和编译选项
ROOT_CFLAGS = $(shell root-config --cflags)
ROOT_LIBS = $(shell root-config --libs)

# 特定的库
SPECTRUM_LIB = -lSpectrum

# 最终的编译和链接选项
FINAL_CXXFLAGS = $(CXXFLAGS) $(ROOT_CFLAGS)
FINAL_LIBS = $(ROOT_LIBS)

# 定义要生成的可执行文件
EXECUTABLES = SSD_calibration apply_calibration

# --- 自动从 config.h 读取输出文件和目录路径 ---
# 使用 sed 从配置文件中提取双引号之间的路径
CONFIG_H = config.h
PARAM_FILE = $(shell grep 'calibration_param_file' $(CONFIG_H) | sed -e 's/.*"\(.*\)"[^"]*$$/\1/')
PLOT_DIR = $(shell grep 'output_plot_dir' $(CONFIG_H) | sed -e 's/.*"\(.*\)"[^"]*$$/\1/')
FINAL_ROOT_FILE = $(shell grep 'final_calibrated_root_file' $(CONFIG_H) | sed -e 's/.*"\(.*\)"[^"]*$$/\1/')


# --- 规则 ---

# "all" 是默认目标，会编译所有程序
all: $(EXECUTABLES)

# 编译 SSD_calibration 的规则
# 新增了 config.h 作为依赖项
SSD_calibration: SSD_calibration.cpp config.h
	$(CXX) $(FINAL_CXXFLAGS) $< -o $@ $(FINAL_LIBS) $(SPECTRUM_LIB)

# 编译 apply_calibration 的规则
# 新增了 config.h 作为依赖项
apply_calibration: apply_calibrationv1.cpp config.h
	$(CXX) $(FINAL_CXXFLAGS) $< -o $@ $(FINAL_LIBS)

# "run" 规则：先编译，然后按顺序执行
run: all
	@echo "--- Running Calibration Step 1: Generating parameters and plots... ---"
	./SSD_calibration
	@echo "\n--- Running Calibration Step 2: Applying parameters to create new ROOT file... ---"
	./apply_calibration
	@echo "\n--- Workflow finished! ---"

# "clean" 规则：删除所有生成的文件
# 使用 -f (force) 和 -r (recursive) 忽略不存在的文件/目录的错误，并递归删除目录
.PHONY: clean all run
clean:
	@echo "Cleaning up compiled executables and object files..."
	rm -f $(EXECUTABLES) *.o
	@echo "Cleaning up generated data files and plots defined in $(CONFIG_H)..."
	rm -f $(PARAM_FILE)
	rm -f $(FINAL_ROOT_FILE)
	rm -rf $(PLOT_DIR)
	@echo "Cleanup complete."