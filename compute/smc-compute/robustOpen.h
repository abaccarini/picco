#pragma once

#include "Operation.h"
#include "poly.hpp"
#include "modMath.hpp"


void RobustOpen(mpz_t result, mpz_t var, bool error_flag, int threadID, NodeNetwork nodeNet, SecretShare *ss) ;
