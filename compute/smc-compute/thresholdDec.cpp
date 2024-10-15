

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
#include "robustOpen.h"

#if __SHAMIR__

using std::cout;
using std::endl;
using std::string;
using std::vector;

/*
contains the thresholdDecryption functionality from Noah's Ark (Dahl et al., 2023)
Parameters:
Q - ciphertext modulus (128-bit large prime)
p - plaintext modulus (2^2)
L - LWE dimension (1024)
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

error_flag determines if we introduce the maximium number of errors into the computation during robustOpening
*/
void SMC_Utils::thresholdDecryption(priv_int *s_prime, mpz_t *a_prime, mpz_t b_prime, int L, int threadID, double &offline_time, double &online_time, bool error_flag) {
    struct timeval tv1;
    struct timeval tv2;
    int _xval = 0;

    gettimeofday(&tv1, NULL);

    mpz_t field;
    mpz_init(field);
    ss->getFieldSize(field);

    // noise term
    mpz_t *E = (mpz_t *)malloc(sizeof(mpz_t) * 2);
    for (int i = 0; i < 2; i++) {
        mpz_init(E[i]);
    }

    mpz_t E_term, v, result, Delta_prime, offset, temp1, temp2, const2;
    mpz_init(E_term);
    mpz_init(v);
    mpz_init(result);
    mpz_init(Delta_prime);
    mpz_init(offset);
    mpz_init(temp1);
    mpz_init(temp2);
    mpz_init_set_ui(const2, 2);

    int log_Bd = 70;
    int pow = 47;
    uint p = 4;
    mpz_fdiv_q_ui(Delta_prime, field, p);

    mpz_pow_ui(offset, const2, pow + log_Bd - 1); // computing 2^(pow + log_Bd - 1)

    // N.A. needs two PRSS's in the range [-2^(pow-1)*Bd , ... , 2^(pow-1)*Bd - 1]
    // Bd \approx 2^70
    // pow = 47
    // => modulus = 2*2^pow*B2 = 2^(pow+1)*Bd
    // subtract generated value by 2^(pow-1)*Bd to obtain value in correct range

    gettimeofday(&tv1, NULL);

    PRSS(pow + log_Bd, 2, E, threadID, ss);
    ss->modSub(E, E, offset, 2);

    // computing [E] <- PRSS.next + PRSS.next
    ss->modAdd(E_term, E_term, E[0]);
    ss->modAdd(E_term, E_term, E[1]);

    gettimeofday(&tv2, NULL);
    offline_time += time_diff(&tv1, &tv2) * 1000.0;

    gettimeofday(&tv1, NULL);
    // computing
    // [v] <- b_prime - a_prime (dot) [s_prime] + [E]
    ss->modDotPub(v, s_prime, a_prime, L);
    ss->modSub(v, b_prime, v);
    ss->modAdd(v, v, E_term);

    RobustOpen(result, v, error_flag, -1, net, ss);
    gettimeofday(&tv2, NULL);

    online_time += time_diff(&tv1, &tv2) * 1000.0;

    for (int i = 0; i < 2; i++)
        mpz_clear(E[i]);
    free(E);

    mpz_clear(E_term);
    mpz_clear(v);
    mpz_clear(result);
    mpz_clear(Delta_prime);
    mpz_clear(offset);
    mpz_clear(temp1);
    mpz_clear(temp2);
    mpz_clear(const2);
}

#endif
