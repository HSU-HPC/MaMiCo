#! /usr/bin/env bash

# Check if source code needs to be formatted using clang-format
# and if CMakeList.txt needs to be formatted with cmake_format

SOURCE_CODE_EXTS=(c h hpp cpp cxx cpph)
CLANG_FORMAT_TARGET=mamico-clangformat

set -e

cd "$(dirname "$0")"
cd ..
CWD=$(pwd)

# Make an ephemeral copy of the working copy
COPY=$(mktemp --directory)
trap 'rm -rf "$COPY"' EXIT
cp -r "$CWD"/* "$COPY"/
cp "$CWD"/.clang-format "$COPY"/.clang-format

# Format files in ephemeral copy
mkdir -p "$COPY/build"
cd "$_"
rm -rf ./*
cmake .. &>/dev/null
make  -j4 "$CLANG_FORMAT_TARGET" &>/dev/null

# Compare source files in working copy and formatted ephemeral copy
find_args=(-type f)
for ext in "${SOURCE_CODE_EXTS[@]}"; do
    find_args+=(-name "*.${ext}" -o)
done
unset 'find_args[-1]'
files=$(find "$CWD" "${find_args[@]}")
for file in $files; do
    relative="${file#"$CWD"/}"
    if ! diff -q "$CWD/$relative" "$COPY/$relative" > /dev/null; then
        echo '**Please format your code with `make '$CLANG_FORMAT_TARGET'`!**' >&2
        exit 1
    fi
done

echo "Code already formatted!"
