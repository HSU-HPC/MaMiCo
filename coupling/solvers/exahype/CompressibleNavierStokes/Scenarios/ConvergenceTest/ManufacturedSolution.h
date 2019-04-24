#ifndef NAVIERSTOKES_MANUFACTUREDSOLUTION__SOURCE__H
#define NAVIERSTOKES_MANUFACTUREDSOLUTION__SOURCE__H

void evaluateSource(double R, double gamma, double kappa, double mu, double t,
                    double x, double y, double *out_3914538540166432036);
void evaluateQ(double gamma, double mu, double t, double x, double y,
               double *out_7076894245671968112);
void evaluateGradQ(double gamma, double mu, double t, double x, double y,
                   double *out_1045687728225358982);

#endif
