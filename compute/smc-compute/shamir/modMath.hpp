#pragma once

#include <gmp.h>

class modMath {
public:
    modMath(mpz_t _PRIME);
    ~modMath() {
        mpz_clear(PRIME);
    };
    modMath(modMath &m) {
        mpz_init(PRIME);
        mpz_set(PRIME, m.PRIME);
    };

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

    // Modular Square root
    void modSqrt(mpz_t, mpz_t);
    void modSqrt(mpz_t *, mpz_t *, int);

    void modSum(mpz_t, mpz_t *, int);

    void copy(mpz_t *src, mpz_t *des, int size);
    void mod(mpz_t *result, mpz_t *a, mpz_t *m, int size);
    void mod(mpz_t *result, mpz_t *a, mpz_t m, int size);

    mpz_t PRIME;
};
