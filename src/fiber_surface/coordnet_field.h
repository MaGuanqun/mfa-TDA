#pragma once
#include "field.h"

namespace fiber_surface {
// Each TorchScript module maps one normalized (t,y,x) point to a scalar.
Field coordnet_field(const std::string& first, const std::string& second, const Box& training_box);
} // namespace fiber_surface
