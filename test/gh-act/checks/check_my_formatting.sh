#! /usr/bin/env bash

# Check if source code needs to be formatted using clang-format
# and if CMakeList.txt needs to be formatted with cmake_format

SOURCE_CODE_EXTS=(c h hpp cpp cxx cpph)
CLANG_FORMAT_TARGET=mamico-clangformat

set -e

cd "$(dirname "$0")"
cd ../../..
CWD=$(pwd)

# Format Python scripts
diff="$(black . --diff 2>/dev/null)"
if [ -n "$diff" ]; then
    echo ':eyes: Please format your Python scripts with `black .`!'
else
    echo ":rocket: Python ccripts already formatted!"
fi

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
    if [[ "$relative" == build/* ]]; then
        continue # Ignore generated files
    fi
    if ! diff -q "$CWD/$relative" "$COPY/$relative" 1>&2; then
        echo ':eyes: Please format your code with `make '$CLANG_FORMAT_TARGET'`!'
        exit 1
    fi
done

echo ":rocket: Code already formatted!"
