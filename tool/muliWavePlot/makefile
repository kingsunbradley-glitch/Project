ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --glibs)

SHOWWAVE_DIR = ../ShowWave_Gui_202511
TARGET = save_wave
OBJECTS = main.o Data_Map.o CAEN_DGTZ.o

CXX = g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++17 -I. -I$(SHOWWAVE_DIR) $(ROOTCFLAGS)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(ROOTLIBS)

main.o: main.cpp $(SHOWWAVE_DIR)/Data_Map.h $(SHOWWAVE_DIR)/CAEN_DGTZ.h $(SHOWWAVE_DIR)/setup.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

Data_Map.o: $(SHOWWAVE_DIR)/Data_Map.cpp $(SHOWWAVE_DIR)/Data_Map.h $(SHOWWAVE_DIR)/setup.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

CAEN_DGTZ.o: $(SHOWWAVE_DIR)/CAEN_DGTZ.cpp $(SHOWWAVE_DIR)/CAEN_DGTZ.h $(SHOWWAVE_DIR)/setup.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)
