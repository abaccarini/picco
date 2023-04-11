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

#ifndef SECRETSHARE_H_
#define SECRETSHARE_H_

#include "SecretShare.h"
#include "stdint.h"
#include <cstdlib>
#include <gmp.h>
#include <iostream>
#include <math.h>
#include <vector>

class SecretShare {

public:
    SecretShare(int, int, mpz_t);

    int getPeers();
    int getThreshold();
    void getFieldSize(mpz_t);

    // Create n shares of a secret or multiple secrets
    void getShares(mpz_t *, mpz_t);
    void getShares(mpz_t **, mpz_t *, int);

    void computeLagrangeWeight();
    void computeSharingMatrix();

    // Reconstruct a secret from n shares
    void reconstructSecret(mpz_t, mpz_t *, bool);
    void reconstructSecret(mpz_t *, mpz_t **, int, bool);

    // Modular Multiplication
    void modMul(mpz_t, mpz_t, mpz_t);
    void modMul(mpz_t *, mpz_t *, mpz_t *, int);
    void modMul(mpz_t, mpz_t, long);
    void modMul(mpz_t *, mpz_t *, long, int);
    void modMul(mpz_t *, mpz_t *, mpz_t, int);

    // Modular Addition
    void modAdd(mpz_t, mpz_t, mpz_t);
    void modAdd(mpz_t *, mpz_t *, mpz_t *, int);
    void modAdd(mpz_t, mpz_t, long);
    void modAdd(mpz_t *, mpz_t *, long, int);
    void modAdd(mpz_t *, mpz_t *, mpz_t, int);

    // Modular Subtraction
    void modSub(mpz_t, mpz_t, mpz_t);
    void modSub(mpz_t *, mpz_t *, mpz_t *, int);
    void modSub(mpz_t, mpz_t, long);
    void modSub(mpz_t, long, mpz_t);
    void modSub(mpz_t *, mpz_t *, long, int);
    void modSub(mpz_t *, long, mpz_t *, int);
    void modSub(mpz_t *, mpz_t *, mpz_t, int);
    void modSub(mpz_t *, mpz_t, mpz_t *, int);

    // Modular Exponentiation
    void modPow(mpz_t, mpz_t, mpz_t);
    void modPow(mpz_t *, mpz_t *, mpz_t *, int);
    void modPow(mpz_t, mpz_t, long);
    void modPow(mpz_t *, mpz_t *, long, int);

    // Modular Inverse
    void modInv(mpz_t, mpz_t);
    void modInv(mpz_t *, mpz_t *, int);

    // Modular Square root
    void modSqrt(mpz_t, mpz_t);
    void modSqrt(mpz_t *, mpz_t *, int);

    // Miscellaneous Functions
    void modSum(mpz_t, mpz_t *, int);
    void copy(mpz_t *src, mpz_t *des, int size);
    void mod(mpz_t *result, mpz_t *a, mpz_t *m, int size);
    void mod(mpz_t *result, mpz_t *a, mpz_t m, int size);

    // computation for 3P multiplication
    void getShares2(mpz_t *temp, mpz_t *rand, mpz_t **data, int size);
    void Seed(unsigned char *key_0, unsigned char *key_1);
    void checkSeed();
    void getCoef(int id);

private:
    mpz_t fieldSize;
    int threshold;
    int peers;
    std::vector<long> coefficients;
    mpz_t *lagrangeWeight;
    mpz_t **sharingMatrix;
    int bits;
    gmp_randstate_t rstate;
    gmp_randstate_t rstate_0;
    gmp_randstate_t rstate_1;
    int seeded;
    int myid;
    int id_p1;
    int id_m1;
    mpz_t id_p1_inv;
};
#endif
