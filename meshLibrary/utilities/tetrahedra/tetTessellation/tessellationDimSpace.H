#ifndef tessellationDimSpace_H
#define tessellationDimSpace_H

#include "direction.H"
#include "scalar.H"

using namespace Foam;

static const scalar VSM = 1e-10;
static const scalar VSM2 = VSM * VSM;

static const direction DIM = 3;
static const direction DIM1 = DIM + 1;
static const direction NV = DIM + 2;

#endif
