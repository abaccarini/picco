/*
   PICCO: A General Purpose Compiler for Private Distributed Computation
   ** Copyright (C) from 2013 PICCO Team
   ** Department of Computer Science and Engineering, University of Notre Dame
   ** Department of Computer Science and Engineering, University of Buffalo (SUNY)

   PICCO is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PICCO is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with PICCO. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef EQZ_SHAMIR_H_
#define EQZ_SHAMIR_H_

#include "Open.h"
#include "Operation.h"
#include "PrefixMultiplication.h"
#include "Random.h"

class EQZ : public Operation {

public:
    EQZ(NodeNetwork nodeNet, std::map<std::string, std::vector<int>> poly, int nodeID, SecretShare *s);
    virtual ~EQZ();
    void doOperation(mpz_t *shares, mpz_t *result, int K, int size, int threaID);

    void doOperation_EQZ(mpz_t *result, mpz_t *a, mpz_t *b, int alen, int blen, int resultlen, int size, int threadID);
    void doOperation_EQZ(mpz_t result, mpz_t a, mpz_t b, int alen, int blen, int resultlen, int threadID);
    void doOperation_EQZ(mpz_t result, mpz_t a, int b, int alen, int blen, int resultlen, int threadID);
private:
    PrefixMultiplication *PreMul;
    Random *Rand;
    mpz_t coef[9];                                      

};

#endif /* EQZ_SHAMIR_H_ */