#pragma once

/* 2024
 * @author Piet Jarmatz
 */

#define DEFINE_DECIMAL_FP6_LIMITS() \
    constexpr double maxFP6 = 1e6; \
    constexpr double stepFP6 = (double)(std::numeric_limits<long long>::max()) / maxFP6; \
    constexpr double minFP6 = 1 / stepFP6; \
    (void)minFP6; // Avoid unused variable warning/error if minFP6 is not needed.

namespace tarch {
namespace utils {

template <class T, class T2> bool contains(T container, T2 element) { return std::find(container.begin(), container.end(), element) != container.end(); }

} // namespace utils
} // namespace tarch
