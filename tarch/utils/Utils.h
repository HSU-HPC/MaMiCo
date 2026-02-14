#pragma once

/* 2024
 * @author Piet Jarmatz
 */

#define DEFINE_DECIMAL_FP_LIMITS(digits) \
    constexpr double maxFP##digits = std::pow(10,digits); \
    constexpr double stepFP##digits = (double)(std::numeric_limits<long long>::max()) / maxFP##digits; \
    constexpr double minFP##digits = 1 / stepFP##digits; \
    (void)minFP##digits; // Avoid unused variable warning/error if minFPx is not needed.

namespace tarch {
namespace utils {

template <class T, class T2> bool contains(T container, T2 element) { return std::find(container.begin(), container.end(), element) != container.end(); }

} // namespace utils
} // namespace tarch
