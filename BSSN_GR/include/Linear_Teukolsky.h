#pragma once

#include "grDef.h"
#include "parameters.h"

namespace bssn {
void LinearTeuk(const double xx, const double yy, const double zz,
                const double t, double *var, bool varsAreGrid = true);
}
