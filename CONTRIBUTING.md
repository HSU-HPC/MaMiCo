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

## External Contribution Policy

These general guidelines should be followed for higher chances of review and merging:
* Make sure the tests pass, and the new code has adequate test coverage
* Use Doxygen annotations. Usage etc. that does not belong in the source code should be covered in the pull request.
* Follow general best practices (descriptive commit messages, multiple smaller commits...)
* Narrow the scope of your pull requests (no kitchen sink approach)
* Code should already be formatted (see above)
* Pull requests should not touch fundamentals (build process, folder structure, core feature, C++ standard etc.)

## AI Usage Policy

* AI contributions must be marked explicitly as AI everywhere, such as in commit messages, source code and in GitHub pull requests and issues.
* In general, we cannot guarantee that AI generated contributions will be considered for review.
