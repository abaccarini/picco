/*
   PICCO: A General Purpose Compiler for Private Distributed Computation
   ** Copyright (C) from 2013 PICCO Team
   ** Department of Computer Science and Engineering, University of Notre Dame
   ** Department of Computer Science and Engineering, University at Buffalo (SUNY)

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

#include "../../../common/shared.h"
#include "../bit_utils.hpp"
#include "ShamirUtil.h"
#include "stdint.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <functional>
#include <gmp.h>
#include <iostream>
#include <map>
#include <math.h>
#include <numeric>
#include <openssl/rand.h>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

// #define KEYSIZE 16

#define COEFF_BOUND 10 // specifies the number of sets of coefs we generate for symmetric function evaluation
#define COEFF_OFFSET 1 // specifies what m value we start with (m = 1), DO NOT CHANGE !!!!!!!!!
#define BASE_10 10     // number base for mpz_t

using std::vector;

class SecretShare {

public:
    SecretShare(unsigned int, unsigned int, mpz_t, unsigned int, unsigned int, unsigned char *[KEYSIZE], std::map<std::string, std::vector<int>>);

    unsigned int getPeers();
    unsigned int getThreshold();
    void getFieldSize(mpz_t);

    unsigned int *getSendToIDs();
    unsigned int *getRecvFromIDs();

    void initCoef();

    void print_poly();

    // Create n shares of a secret or multiple secrets
    void getShares(mpz_t *, mpz_t);
    void getShares(mpz_t **, mpz_t *, int);

    void computeLagrangeWeights();
    void computeSharingMatrix();

    void computeLagrangePolys();

    // Reconstruct a secret from n shares
    void reconstructSecret(mpz_t, mpz_t *);
    void reconstructSecret(mpz_t *, mpz_t **, int);
    // Reconstruct a secret from the minimum (threshold+1) number of shares
    void reconstructSecretFromMin(mpz_t *, mpz_t **, unsigned int);
    void reconstructSecretMult(mpz_t *result, mpz_t **y, int size);
    void reconstructSecretFromMin_test(mpz_t *, mpz_t **, unsigned int);

    // Evaluate a polynomial represented by threshold+1 shares on another threshold+1 points
    void getSharesMul(mpz_t **, mpz_t **, unsigned int);

    void modDotPub(mpz_t result, mpz_t *x, mpz_t *pub, int size);
    void modDotPub(mpz_t result, mpz_t *x, int *y, int size);

    // Modular Multiplication
    void modMul(mpz_t, mpz_t, mpz_t);
    void modMul(mpz_t *, mpz_t *, mpz_t *, int);
    void modMul(mpz_t, mpz_t, long);
    void modMul(mpz_t *, mpz_t *, long, int);
    void modMul(mpz_t *, mpz_t *, mpz_t, int);
    void modMul(mpz_t *result, mpz_t *x, int *y, int size);

    // Modular Addition
    void modAdd(mpz_t, mpz_t, mpz_t);
    void modAdd(mpz_t *, mpz_t *, mpz_t *, int);
    void modAdd(mpz_t, mpz_t, long);
    void modAdd(mpz_t *, mpz_t *, long, int);
    void modAdd(mpz_t *, mpz_t *, mpz_t, int);
    void modAdd(mpz_t *, mpz_t *, long *, int);
    void modAdd(mpz_t *, mpz_t *, int *, int);

    // Modular Subtraction
    void modSub(mpz_t, mpz_t, mpz_t);
    void modSub(mpz_t *, mpz_t *, mpz_t *, int);
    void modSub(mpz_t, mpz_t, long);
    void modSub(mpz_t, long, mpz_t);
    void modSub(mpz_t *, mpz_t *, long, int);
    void modSub(mpz_t *, long, mpz_t *, int);
    void modSub(mpz_t *, mpz_t *, mpz_t, int);
    void modSub(mpz_t *, mpz_t, mpz_t *, int);
    void modSub(mpz_t *result, mpz_t *x, int *y, int size);
    void modSub(mpz_t *result, int *x, mpz_t *y, int size);

    // Modular Exponentiation
    void modPow(mpz_t, mpz_t, mpz_t);
    void modPow(mpz_t *, mpz_t *, mpz_t *, int);
    void modPow(mpz_t, mpz_t, long);
    void modPow(mpz_t *, mpz_t *, long, int);

    void modPow2(mpz_t, int);
    void modPow2(mpz_t, mpz_t);
    void modPow2(mpz_t *result, int *exponent, int size);
    void modPow2(mpz_t *result, mpz_t *exponent, int size);

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

    // computation for 3P multiplication, to be removed
    void getShares2(mpz_t *temp, mpz_t *rand, mpz_t **data, int size);
    void Seed(unsigned char *key_0, unsigned char *key_1);
    void getCoef(int id);

    void PRG(mpz_t **output, uint size, uint start_ind);
    void PRG_thread(mpz_t **output, uint size, uint start_ind, int threadID);

    int computePolynomials(std::vector<int> polys, int point);

    // uint computeMultPolySize(uint A_size, uint B_size);
    // int computePolyLen(mpz_t *A, int A_len);
    // void polyMult(mpz_t *result, mpz_t *A, mpz_t *B, uint A_size, uint B_size);
    // void polyDivRemainder(mpz_t *quotient, mpz_t *remainder, mpz_t *A, mpz_t *B, uint A_size, uint B_size, uint &Q_size, uint &R_size);
    // void interpolate(mpz_t *result, mpz_t *y);
    // void gaoDecoding(mpz_t *result, mpz_t *shares);
    void randInit(unsigned char *keys[KEYSIZE]);
    void randInit_thread(int threadID);
    void randInit_thread_mult(int threadID);

    void generateRandValue(int bits, int size, mpz_t *results);
    void generateRandValue(int bits, int size, mpz_t *results, int threadID);
    void generateRandValue(mpz_t mod, int size, mpz_t *results);
    void generateRandValue(mpz_t mod, int size, mpz_t *results, int threadID);

    void PRZS(mpz_t mod, int size, mpz_t *results);

    int getCoefIndex(int k);

    std::vector<std::string> splitfunc(const char *str, const char *delim);
    std::vector<std::string> split(const std::string s, const std::string delimiter, int expected_size = 0);

    bool is_int(const std::string &str);
    bool is_float(const std::string &str);

    void ss_input(int id, int *var, std::string type, std::ifstream *inputStreams);
    void ss_input(int id, mpz_t *var, std::string type, std::ifstream *inputStreams);
    void ss_input(int id, float *var, std::string type, std::ifstream *inputStreams);
    void ss_input(int id, mpz_t **var, std::string type, std::ifstream *inputStreams);
    void ss_input(int id, mpz_t *var, int size, std::string type, std::ifstream *inputStreams);
    void ss_input(int id, int *var, int size, std::string type, std::ifstream *inputStreams);
    void ss_input(int id, mpz_t **var, int size, std::string type, std::ifstream *inputStreams);
    void ss_input(int id, float *var, int size, std::string type, std::ifstream *inputStreams);

    void ss_input(mpz_t *var, std::string type);
    void ss_input(mpz_t **var, std::string type);
    void ss_input(mpz_t *var, int size, std::string type);
    void ss_input(mpz_t **var, int size, std::string type);

    void ss_output(int id, int *var, std::string type, std::ofstream *outputStreams);
    void ss_output(int id, mpz_t *var, std::string type, std::ofstream *outputStreams);
    void ss_output(int id, float *var, std::string type, std::ofstream *outputStreams);
    void ss_output(int id, mpz_t **var, std::string type, std::ofstream *outputStreams);
    void ss_output(int id, mpz_t *var, int size, std::string type, std::ofstream *outputStreams);
    void ss_output(int id, int *var, int size, std::string type, std::ofstream *outputStreams);
    void ss_output(int id, mpz_t **var, int size, std::string type, std::ofstream *outputStreams);
    void ss_output(int id, float *var, int size, std::string type, std::ofstream *outputStreams);

    void ss_single_convert_to_private_float(float a, mpz_t **priv_a, int len_sig, int len_exp);

    void ss_process_operands(mpz_t **a1, mpz_t **b1, int alen_sig, int alen_exp, int blen_sig, int blen_exp, int *len_sig, int *len_exp, int size);
    mpz_t **coef; // public because eqz uses it directly

private:
    std::map<std::string, std::vector<int>> polynomials; // public for easier access in Random, but polynomials are only accessed inside of generateRandomValue?

    mpz_t fieldSize;
    unsigned int threshold;
    unsigned int peers;
    unsigned int myID;

    unsigned int numThreads;

    std::vector<long> coefficients;
    // coeffiicents for polynomial reconstruction on point 0 from all shares
    mpz_t *lagrangeWeightsAll;
    // coefficients for polynomial reconstruction on point 0 from threshold+1 shares at points with indices in recvFromIDs
    mpz_t *lagrangeWeightsThreshold;
    // coefficients for polynonial evaluation on threshold points at indices stored in sendToIDs, where the polynomial is stored as threshold+1 values at indices in recvFromIDs and point 0
    // size is threshold*(threshold+1)
    mpz_t **lagrangeWeightsMult;
    mpz_t **sharingMatrix;
    //    int bits;
    gmp_randstate_t rstate;

    mpz_t **lagrangePolysAll; // used for robust reconstruction

    // peers to whom a share or shares will be sent, numbered consequently
    // from myID (myID+1, ..., myID+t)
    // void checkSeed();
    unsigned int *sendToIDs;
    // peers to receive shares from or generate via PRGs, numbered from myID
    // in the decreasing order (myID-t, ..., myID-1) ***this is ultimately INCREASING order
    // e.g. for id = 3, recvFromIDs[0] = 4, recvFromIDs[1] = 5
    unsigned int *recvFromIDs;

    uint *multIndices;

    unsigned char **mult_keys;
    // additional data structures for multiplication
    static gmp_randstate_t *rstatesMult;

    static int *rand_isFirst_thread_mult;
    gmp_randstate_t rstate_mine;
    static gmp_randstate_t **rstates_thread_mult;

    // for 3-party multiplication
    gmp_randstate_t rstate_0;
    gmp_randstate_t rstate_1;
    int seeded;
    int id_p1;
    int id_m1;
    mpz_t id_p1_inv;

    // from Random.cpp
    // int rand_isInitialized;
    static int *rand_isFirst_thread;
    static gmp_randstate_t *rstates;
    static gmp_randstate_t **rstates_thread;
    static pthread_mutex_t mutex;

    mpz_t *poly_evaluation; // stores polynomial evaluation for random generation, only is filled once
};

vector<int> generateCoef(int m, uint &inv_term);
vector<int> multiply_poly(vector<int> A, vector<int> B);

template <typename T>
int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

uint gcd(uint a, uint b);

uint lcm(uint a, uint b);

// void ss_clear(mpz_t &x);
// void ss_set_str(mpz_t x, const char *str, int base);

// char *ss_get_str(char *str, int base, const mpz_t op);

// void ss_free_arr(mpz_t *op, int size);
// void ss_init_set_si(mpz_t &x, int x_val);

#endif
