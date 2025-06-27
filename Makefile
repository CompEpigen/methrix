# Simple Makefile for static ModKit merge tool
# Usage: module load HTSlib/1.20-GCC-14.1.0 && make

# Compiler and flags
CC = gcc
CFLAGS = -O3 -Wall -std=c99
TARGET = modkit_merge_static
SOURCE = modkit_merge_improved.c

# Check if HTSlib module is loaded
ifndef EBROOTHTSLIB
$(error HTSlib module not loaded. Run: module load HTSlib/1.20-GCC-14.1.0)
endif

# HTSlib paths from module
HTSLIB_INCLUDE = -I$(EBROOTHTSLIB)/include
HTSLIB_LIB = -L$(EBROOTHTSLIB)/lib

# Static linking strategy: HTSlib + compression libs static, system libs dynamic
STATIC_LIBS = -Wl,-Bstatic -lhts -lz -lbz2 -llzma
DYNAMIC_LIBS = -Wl,-Bdynamic -lm -lpthread -ldl

.PHONY: all clean test check-deps static

all: $(TARGET)

$(TARGET): $(SOURCE)
@echo "Building static binary with selective linking..."
@echo "HTSlib: $(EBROOTHTSLIB)"
$(CC) $(CFLAGS) $(HTSLIB_INCLUDE) -o $(TARGET) $(SOURCE) \
$(HTSLIB_LIB) $(STATIC_LIBS) $(DYNAMIC_LIBS)
@echo "Static binary created: $(TARGET)"

static: $(TARGET)
@echo "Checking dependencies..."
@ldd $(TARGET) || echo "Fully static binary (no dependencies)"
@echo "Binary size: $$(du -h $(TARGET) | cut -f1)"

test: $(TARGET)
@echo "Testing static binary..."
@./$(TARGET) --help > /dev/null && echo "✓ Binary works!" || echo "✗ Binary failed"

check-deps: $(TARGET)
@echo "=== Dependency Analysis ==="
@echo "Dynamic dependencies:"
@ldd $(TARGET)
@echo ""
@echo "Binary size: $$(du -h $(TARGET) | cut -f1)"
@echo "Expected: Only system libs (libc, libm, libpthread) should be dynamic"

clean:
  rm -f $(TARGET)

# Help target
help:
  @echo "ModKit Merge Static Build"
@echo ""
@echo "Prerequisites:"
@echo "  module load HTSlib/1.20-GCC-14.1.0"
@echo ""
@echo "Targets:"
@echo "  all        - Build static binary (default)"
@echo "  static     - Build and check dependencies"
@echo "  test       - Test the binary"
@echo "  check-deps - Analyze binary dependencies"
@echo "  clean      - Remove binary"
@echo "  help       - Show this help"
@echo ""
@echo "Usage:"
@echo "  make && make check-deps"