# Contributing Guidelines

## Branch Naming

Use the suffixes `-WIP` or `_WIP` if you **don't** want tests to run on your branch before opening a pull request.

Consider using a sensible branch prefix like `fix/`, `feature/` or `enhance/`.

## Code Formatting

Format your code, CMakeList.txt, and Python scripts using the following commands:

|Files|Tool|Version|Command|
|:---|:---|---|:---|
|C/C++ source|[clang-format](https://releases.llvm.org/18.1.8/tools/clang/docs/ClangFormat.html)|18.1.8|`make mamico-clangformat`|
|CMakeList.txt|[cmakelang](https://github.com/cheshirekow/cmake_format)|0.6.13|`cmake-format --config-file .cmake-format.py -i CMakeLists.txt`|
|Python scripts|[black](https://github.com/psf/black)|26.5.1|`black .`|
