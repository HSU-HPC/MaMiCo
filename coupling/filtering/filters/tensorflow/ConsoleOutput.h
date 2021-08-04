#pragma once

#include <vector>
#include <iostream>

//prints 1D vector contents to console
template <class A>
std::ostream& operator<< (std::ostream& s, const std::vector<A> vec){
	for(auto entry : vec) s<<entry<<"  ";
	return s;
}

//prints 2D vector contents to console
template <class A>
std::ostream& operator<< (std::ostream& s, const std::vector<std::vector<A>> vec){
	for(auto entry : vec) s<<entry<<"\n"<<std::endl;
	return s;
}
