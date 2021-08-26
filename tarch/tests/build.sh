g++ -c ../tinyxml2/tinyxml2.cpp -o tinyxml2.o
g++ main_readconfig.cpp tinyxml2.o -o readconfig -I../../ -std=c++0x

