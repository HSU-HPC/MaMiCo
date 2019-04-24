#ifndef COMPRESSIBLENAVIERSTOKES_VARIANCEHELPER_H
#define COMPRESSIBLENAVIERSTOKES_VARIANCEHELPER_H

#include <utility>

std::pair<double, double> mergeVariance(double mean0, double mean1, double var0,
                                        double var1, double count0, double count1);

#endif  // COMPRESSIBLENAVIERSTOKES_VARIANCEHELPER_H
