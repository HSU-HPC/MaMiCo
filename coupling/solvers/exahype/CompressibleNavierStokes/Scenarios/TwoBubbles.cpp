#include "TwoBubbles.h"

NavierStokes::TwoBubbles::TwoBubbles() {
  // Robert (1993),
  // https://doi.org/10.1175/1520-0469(1993)050<1865:BCEWAS>2.0.CO;2
  const double A_1 = 0.5;  // temp. difference
  const double a_1 = 150;  // [m]
  const double S_1 = 50;   // [m]
  const double x_1 = 500;  // [m]
  const double z_1 = 300;  // [m]

  // Small, cold bubble
  const double A_2 = -0.15;  // temp. difference
  const double a_2 = 0;      // [m]
  const double S_2 = 50;     // [m]
  const double x_2 = 560;    // [m]
  const double z_2 = 640;    // [m]

  const auto bubbleType = BubbleType::smooth;

  bubbles = std::vector<CloudScenario::Bubble>{
      Bubble(bubbleType, A_1, a_1, S_1, x_1, z_1),
      Bubble(bubbleType, A_2, a_2, S_2, x_2, z_2),
  };
}
