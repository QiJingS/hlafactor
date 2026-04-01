HTSROOT = $(CURDIR)/_install

CXXFLAGS += -I$(HTSROOT)/include
LDFLAGS  += -L$(HTSROOT)/lib
LDLIBS   += -lhts
CXX = g++
CXXFLAGS = -std=c++20 -O2 -I./htslib
LDFLAGS = -L./htslib -lhts

TARGET = ../hlafactor
SRC = hlafactor.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS) -Wl,-rpath,@loader_path/htslib

clean:
	rm -f $(TARGET)

.PHONY: all clean