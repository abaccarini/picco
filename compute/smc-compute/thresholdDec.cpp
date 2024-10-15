

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

#include "SMC_Utils.h"
#include "ShamirUtil.h"

#if __SHAMIR__

using std::cout;
using std::endl;
using std::string;
using std::vector;

/*
contains the thresholdDecryption functionality from Noah's Ark (Dahl et al., 2023)
Parameters:
Q - ciphertext modulus
p - plaintext modulus
L - LWE dimension
three values are entered into the computation:
    s_prime: a secret-shared L-bit secret key, used to compute b_prime
    a_prime: a *public* random L-vector (generated elsewhere)
    b_prime: the computed ciphertext of the form:

        b_prime = a_prime (dot) s_prime + e + \Delta*m

    Delta = floor(Q/p), e - some (small) discrete gaussian noise

Decryption is performed as follows (arithmetic modulo F_Q):
    - Generate two PRSS's in the range [-2^(pow-1)*Bd , ... , 2^(pow-1)*Bd - 1] and sum together to [E]
    - compute [v] <- b_prime - a_prime (dot) [s_prime] + [E]
    - Call RobustOpen([v]) to obtain v = b -  a_prime (dot) s_prime + E
    - Recover m by computing round(v/Delta) (mod p)
*/
void SMC_Utils::thresholdDecryption(priv_int *s_prime, mpz_t *a_prime, mpz_t b_prime, int size, int threadID) {

    // noise term
    mpz_t *E = (mpz_t *)malloc(sizeof(mpz_t) * 2);

    for (int i = 0; i < 2; i++) {
        mpz_init(E[i]);
    }

    mpz_t E_term;
    mpz_init(E_term);

    mpz_t temp1, temp2, const1, const2, constK, constK_m1, constm, const2K_m1, const2K, const2m;
    mpz_init(temp1);
    mpz_init(temp2);
    mpz_init_set_ui(const1, 1);
    mpz_init_set_ui(const2, 2);
    int log_Bd = 70;
    int pow = 47;
    // mpz_init_set_ui(constK, K);
    // mpz_init_set_ui(constK_m1, K - 1);
    // mpz_init_set_ui(constm, m);

    mpz_t offset;
    mpz_init(offset);
    mpz_pow_ui(offset, const2, pow + log_Bd - 1);
    // mpz_t modulus, pow_term, Bd;
    // mpz_init(pow_term);
    // mpz_init(Bd);

    // mpz_init(const2K);
    // mpz_init(const2K_m1);
    // mpz_init(const2m);

    // mpz_pow_ui(pow_term, const2, 47); // 2^47
    // // mpz_pow_ui(pow_term, const2, 47); // 2^47
    // mpz_sub(temp1, pow_term, const1);
    // mpz_mul(temp1, temp1, Bd); // temp1 contains (2^pow -1)*2^70

    // long long nu = nChoosek(ss->getPeers(), ss->getThreshold());
    // mpz_fdiv_q_ui(modulus, temp1, nu); // computes ((2^pow -1)*2^70) / (n Choose t)

    // N.A. needs two PRSS's in the range [-2^(pow-1)*Bd , ... , 2^(pow-1)*Bd - 1]
    // Bd \approx 2^70
    // pow = 47
    // => modulus = 2*2^pow*B2 = 2^(pow+1)*Bd
    // subtract generated value by 2^(pow-1)*Bd to obtain value in correct range
    PRSS(pow + log_Bd, 2, E, threadID, ss);
    ss->modSub(E, E, offset, 2);

    ss->modAdd(E_term, E_term, E[0]);
    ss->modAdd(E_term, E_term, E[1]);



    for (int i = 0; i < 2; i++)
        mpz_clear(E[i]);
    free(E);
}

#endif
