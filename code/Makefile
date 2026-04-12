UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
include Makefile.macos
else
include Makefile.linux
endif