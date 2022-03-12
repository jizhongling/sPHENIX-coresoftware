#include "TrainingHits.h"

TrainingHits::TrainingHits()
{
  Reset();
}

void TrainingHits::Reset()
{
  v_adc.fill(0);
  phi = 0.;
  z = 0.;
  layer = 0;
}
