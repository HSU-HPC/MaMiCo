#pragma once

#include <vector>
#include <iostream>

template <class A>
std::ostream& operator<< (std::ostream& s, const std::vector<A> vec);

template <class A>
std::ostream& operator<< (std::ostream& s, const std::vector<std::vector<A>> vec);

//#include "ConsoleOutput.cpph"