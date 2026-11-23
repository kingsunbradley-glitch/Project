CXX ?= g++
ROOT_CONFIG ?= root-config
ANA_SETUP_DIR ?= /home/evalie/Project/Function/ANA/ANA_SHANS2_202602_v1.07.4/Cal_Auto_All_202412/Cal_Auto_DSSD_Normalize
TARGET ?= Normalize_DSSD

CPPFLAGS += -I"$(ANA_SETUP_DIR)" -DNORM_SETUP_PATH='"$(ANA_SETUP_DIR)/setup.h"'
CXXFLAGS += -O2 -Wall -Wextra $(shell $(ROOT_CONFIG) --cflags)
LDLIBS += $(shell $(ROOT_CONFIG) --libs) -lMinuit2

.PHONY: all clean
all: $(TARGET)

$(TARGET): main.cpp RobustFit.cpp RobustFit.h NormConfig.h $(ANA_SETUP_DIR)/setup.h makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o "$@" main.cpp RobustFit.cpp $(LDLIBS)

clean:
	rm -f "$(TARGET)"
