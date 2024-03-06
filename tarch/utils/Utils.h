#pragma once

/* 2024
 * @author Piet Jarmatz
 */

namespace tarch {
namespace utils {

template<class T, class T2> bool contains(T container, T2 element){
	return std::find(container.begin(), container.end(), element) != container.end();
}

}
}
