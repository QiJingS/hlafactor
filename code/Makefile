CXX ?= g++

APP := hlafactor

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
PLATFORM := $(UNAME_S)-$(UNAME_M)

CURDIR_ABS := $(CURDIR)
ROOT_DIR := $(CURDIR_ABS)/..
HTSLIB_SRC ?= $(CURDIR_ABS)/htslib

PREFIX := $(CURDIR_ABS)/_install/$(PLATFORM)
BUILD_DIR := $(CURDIR_ABS)/build/$(PLATFORM)
DIST_DIR  := $(CURDIR_ABS)/dist

VERSION ?= 0.1.0
PKG_DIR := $(DIST_DIR)/$(APP)-$(VERSION)-$(PLATFORM)

SRCS := \
	hlafactor.cpp

OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.o)

CPPFLAGS += -I$(PREFIX)/include -I$(HTSLIB_SRC)
CXXFLAGS += -std=c++20 -O2 -Wall -Wextra

LDFLAGS += -L$(PREFIX)/lib
LDLIBS   += -lhts -lz -llzma -lm -lpthread

ifeq ($(UNAME_S),Darwin)
	RPATH := -Wl,-rpath,@executable_path/../lib
else
	RPATH := -Wl,-rpath,'$$ORIGIN/../lib'
endif

all: $(BUILD_DIR)/$(APP)

deps: $(PREFIX)/.built

$(PREFIX)/.built:
	@echo "Building htslib for $(PLATFORM) into $(PREFIX)"
	@mkdir -p $(PREFIX)
	cd $(HTSLIB_SRC) && ./configure --prefix="$(PREFIX)" --disable-bz2 && $(MAKE) -j && $(MAKE) install
	@touch $@

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/$(APP): $(OBJS) $(PREFIX)/.built | $(BUILD_DIR)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS) $(RPATH)

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

package: all
	rm -rf $(PKG_DIR)
	mkdir -p $(PKG_DIR)/bin
	mkdir -p $(PKG_DIR)/lib
	mkdir -p $(PKG_DIR)/data_input

	cp $(BUILD_DIR)/$(APP) $(PKG_DIR)/bin/$(APP).bin
	printf '%s\n' \
		'#!/bin/sh' \
		'SCRIPT_DIR="$$(CDPATH= cd -- "$$(dirname -- "$$0")" && pwd)"' \
		'exec "$$SCRIPT_DIR/$(APP).bin" "$$@"' \
		> $(PKG_DIR)/bin/$(APP)
	chmod +x $(PKG_DIR)/bin/$(APP)

	cp -a $(ROOT_DIR)/data_input/. $(PKG_DIR)/data_input/
	cp $(ROOT_DIR)/README.md $(PKG_DIR)/
	cp -a $(PREFIX)/lib/libhts.* $(PKG_DIR)/lib/ 2>/dev/null || true

	tar -czf $(PKG_DIR).tar.gz -C $(DIST_DIR) $(notdir $(PKG_DIR))

clean:
	rm -rf $(BUILD_DIR)

distclean: clean
	rm -rf $(PREFIX) $(DIST_DIR)

.PHONY: all deps package clean distclean