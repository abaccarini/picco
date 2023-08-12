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
#include "PrefixMultiplication.h"

PrefixMultiplication::PrefixMultiplication(NodeNetwork nodeNet, std::map<std::string, std::vector<int>> poly, int nodeID, SecretShare *s) {

    net = nodeNet;
    polynomials = poly;
    id = nodeID;
    ss = s;

    Rand = new Random(nodeNet, poly, nodeID, s);
}

PrefixMultiplication::~PrefixMultiplication() {}

// Source: Catrina and de Hoogh, "Improved Primites for Secure Multiparty Integer Computation," 2010
// Protocol 4.2 page 11
// B[length][size]
void PrefixMultiplication::doOperation(mpz_t **B, mpz_t **result, int length, int size, int threadID) {
    int peers = ss->getPeers();
    mpz_t *R = (mpz_t *)malloc(sizeof(mpz_t) * length * size);
    mpz_t *S = (mpz_t *)malloc(sizeof(mpz_t) * length * size);
    mpz_t *U = (mpz_t *)malloc(sizeof(mpz_t) * length * size);
    mpz_t *V = (mpz_t *)malloc(sizeof(mpz_t) * length * size);
    // mpz_t *results = (mpz_t *)malloc(sizeof(mpz_t) * length * size);

    // these are needed for Step 5 multiplication
    mpz_t **buffer1 = (mpz_t **)malloc(sizeof(mpz_t *) * peers);
    mpz_t **buffer2 = (mpz_t **)malloc(sizeof(mpz_t *) * peers);

    for (int i = 0; i < length * size; i++) {
        mpz_init(R[i]);
        mpz_init(S[i]);
        mpz_init(V[i]);
        mpz_init(U[i]);
        // mpz_init(results[i]);
    }

    for (int i = 0; i < size; i++)
        for (int i = 0; i < peers; i++) {
            buffer1[i] = (mpz_t *)malloc(sizeof(mpz_t) * length * size);
            buffer2[i] = (mpz_t *)malloc(sizeof(mpz_t) * length * size);
            for (int j = 0; j < length * size; j++) {
                mpz_init(buffer1[i][j]);
                mpz_init(buffer2[i][j]);
            }
        }

    mpz_t field;
    mpz_init(field);
    ss->getFieldSize(field);
    Rand->generateRandValue(id, field, length * size, R, threadID); //
    Rand->generateRandValue(id, field, length * size, S, threadID); //

    // step 4, MulPub (can't be replaced with Open)
    ss->modMul(U, R, S, length * size);
    net.broadcastToPeers(U, length * size, buffer1, threadID);
    ss->reconstructSecret(U, buffer1, length * size);
    clearBuffer(buffer1, peers, length * size);

    // step 5, multiplication (not using the mult object?) -- i think its because theres this offset by one of R, whoever originally wrote this's logic
    // (length - 1) * size total multiplications
    for (int i = 0; i < length - 1; i++)
        for (size_t j = 0; j < size; j++) {
            // ss->modMul(V[i * size + j], R[i * size + j + 1], S[i * size + j]);
            // i think this is the correct version?
            ss->modMul(V[j * size + i], R[j * size + i + 1], S[j * size + i]);
        }
    ss->getShares(buffer1, V, (length - 1) * size);
    net.multicastToPeers(buffer1, buffer2, (length - 1) * size, threadID);
    ss->reconstructSecret(V, buffer2, (length - 1) * size);
    clearBuffer(buffer1, peers, length * size);
    // end step 5

    // computing all the inverses of u (used in steps 7 and 8, only need to be done once)
    ss->modInv(U, U, length * size);

    // mpz_set(W[0], R[0]); // not needed since we are using R in place of U
    for (int i = 1; i < length; i++) {
        // step 7
        ss->modMul(R[i], V[i - 1], U[i - 1]);
    }

    // step 8
    ss->modMul(S, S, U, length * size);

    // step 9, MulPub (can't be replaced with Open)
    for (int i = 0; i < length; i++)
        for (int m = 0; m < size; m++)
            // computing m_i <- [a_i] * [w_i]
            ss->modMul(U[i * size + m], B[i][m], R[i * size + m]); // CHECK INDICES FOR U/R/B
    net.broadcastToPeers(U, size * length, buffer1, threadID);
    ss->reconstructSecret(U, buffer1, size * length);
    // end step 9

    for (int i = 0; i < size; i++)
        mpz_set(R[i], U[i]); // reusing R, setting all the m_1's to first (size) elements of R

    ss->copy(B[0], result[0], size); // all the a_0's are stored in b[0]
    for (int i = 1; i < length; i++) {
        for (int j = 0; j < size; j++) {
            // computing prod(m_i, i = 1, j)
            ss->modMul(R[j], R[j], U[i * size + j]);
            // computing [p_j] = [z_j] * prod(m_i, i = 1, j)
            ss->modMul(result[i][j], S[i], R[j]);
        }
    }

    for (int i = 0; i < peers; i++) {
        for (int j = 0; j < length; j++) {
            mpz_clear(buffer1[i][j]);
            mpz_clear(buffer2[i][j]);
        }
        free(buffer1[i]);
        free(buffer2[i]);
    }
    free(buffer1);
    free(buffer2);

    for (int i = 0; i < length; i++) {
        mpz_clear(R[i]);
        mpz_clear(S[i]);
        mpz_clear(V[i]);
        mpz_clear(U[i]);
    }
    free(R);
    free(S);
    free(V);
    free(U);

    // for (int i = 0; i < length * size; i++)
    //     mpz_clear(results[i]);
    // free(results);
}
