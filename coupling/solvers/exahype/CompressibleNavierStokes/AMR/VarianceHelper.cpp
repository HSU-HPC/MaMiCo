#include "VarianceHelper.h"

std::pair<double, double> mergeVariance(double mean0, double mean1, double var0,
                                        double var1, double count0, double count1) {
  // Merge using Chan et al. algorithm
  // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Parallel_algorithm
  if (count0 == 0) {
    return {mean1, var1};
  }
  if (count1 == 0) {
    return {mean0, var0};
  }

  const auto totalCount = count0 + count1;

  const auto delta = mean1 - mean0;
  const auto mA = var0 * (count0 - 1);
  const auto mB = var1 * (count1 - 1);
  const auto mTotal =
      mA + mB + (delta * delta * count0 * count1) / (totalCount);

  const auto mergedMean = (mean0 * count0 + mean1 * count1) / (totalCount);
  const auto mergedVariance = mTotal / (totalCount - 1);
  // Note: MergedVariance uses Bessel's correction!

  return {mergedMean, mergedVariance};
}
