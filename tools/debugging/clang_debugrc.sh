#!/bin/bash

# Check if the script is being sourced
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR! This script must be sourced, not executed."
    echo "Use: source $0"
    exit 1
fi

# See https://clang.llvm.org/docs/AddressSanitizer.html
export CXX=clang++
export CC=clang
export ASAN_OPTIONS="detect_leaks=1:detect_stack_use_after_return=1:check_initialization_order=1:strict_init_order=1:alloc_dealloc_mismatch=1:halt_on_error=1:abort_on_error=1"
export UBSAN_OPTIONS="print_stacktrace=1:halt_on_error=1:abort_on_error=1"
if command -v llvm-symbolizer >/dev/null 2>&1; then
    ASAN_SYMBOLIZER_PATH="$(command -v llvm-symbolizer)"
    export ASAN_SYMBOLIZER_PATH
fi

echo "Clang debug environment loaded."