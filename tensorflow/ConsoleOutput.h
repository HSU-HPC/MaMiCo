
#pragma once

#include <vector>
#include <iostream>

template <class A>
std::ostream& operator<< (std::ostream& s, const std::vector<A> vec){
	for(auto entry : vec) s<<entry<<"  ";
	return s;
}

template <class A>
std::ostream& operator<< (std::ostream& s, const std::vector<std::vector<A>> vec){
	for(auto entry : vec) s<<entry<<"\n"<<std::endl;
	return s;
}
