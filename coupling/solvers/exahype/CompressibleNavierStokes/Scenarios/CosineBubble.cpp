#include "CosineBubble.h"

NavierStokes::CosineBubble::CosineBubble() {
  const auto tempDiff = 0.5;
  const auto x = 500;     // [m]
  const auto z = 350;     // [m]
  const auto size = 250;  // [m]
  const auto decay = -1;  // doesn't matter for cosine bubble

  const auto bubbleType = BubbleType::cosine;

  bubbles = std::vector<CloudScenario::Bubble>{
      Bubble(bubbleType, tempDiff, size, decay, x, z),
  };
}
