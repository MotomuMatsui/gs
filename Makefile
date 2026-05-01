# ===========================================================================
# Graph Splitting Method v2.5 — top-level Makefile (platform dispatcher)
#
# Select the Makefile for your platform:
#
#   Linux  (GCC, optional NVIDIA CUDA):
#     make -f Makefile.linux           # CPU-only  → gs2
#     make -f Makefile.linux gpu       # GPU build → gs2_gpu
#     make -f Makefile.linux clean_all
#
#   macOS  (Apple clang + Homebrew libomp, CPU-only):
#     brew install gcc libomp cmake    # one-time setup
#     make -f Makefile.macos           # CPU-only  → gs2
#     make -f Makefile.macos clean_all
#
# This file itself just forwards 'make' to the correct platform Makefile.
# ===========================================================================

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
  PLATFORM_MK := Makefile.macos
else
  PLATFORM_MK := Makefile.linux
endif

# Forward every target to the platform-specific Makefile
%:
	$(MAKE) -f $(PLATFORM_MK) $@

.PHONY: all
all:
	$(MAKE) -f $(PLATFORM_MK) all
