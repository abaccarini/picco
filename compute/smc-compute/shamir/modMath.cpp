
#include "modMath.hpp"
#include <cstdlib>
#include <iostream>
modMath::modMath(mpz_t _PRIME) {
    mpz_init(PRIME);
    mpz_set(PRIME, _PRIME);
}

void modMath::modMul(mpz_t result, mpz_t x, mpz_t y) {
    mpz_mul(result, x, y);
    mpz_mod(result, result, PRIME);
}

void modMath::modMul(mpz_t *result, mpz_t *x, mpz_t *y, int size) {
    for (int i = 0; i < size; i++)
        modMul(result[i], x[i], y[i]);
}

void modMath::modMul(mpz_t result, mpz_t x, long y) {
    mpz_mul_si(result, x, y);
    mpz_mod(result, result, PRIME);
}

void modMath::modMul(mpz_t *result, mpz_t *x, long y, int size) {
    for (int i = 0; i < size; ++i) {
        modMul(result[i], x[i], y);
    }
}

void modMath::modMul(mpz_t *result, mpz_t *x, mpz_t y, int size) {
    for (int i = 0; i < size; ++i) {
        modMul(result[i], x[i], y);
    }
}

void modMath::modMul(mpz_t *result, mpz_t *x, int *y, int size) {
    mpz_t *ytmp = (mpz_t *)malloc(sizeof(mpz_t) * size);
    for (int i = 0; i < size; i++) {
        mpz_init_set_si(ytmp[i], y[i]);
        modMul(result[i], ytmp[i], x[i]);
    }
    for (int i = 0; i < size; i++)
        mpz_clear(ytmp[i]);
}

void modMath::modAdd(mpz_t result, mpz_t x, mpz_t y) {
    mpz_add(result, x, y);
    mpz_mod(result, result, PRIME);
}

void modMath::modAdd(mpz_t *result, mpz_t *x, mpz_t *y, int size) {
    for (int i = 0; i < size; i++)
        modAdd(result[i], x[i], y[i]);
}

void modMath::modAdd(mpz_t result, mpz_t x, long y) {
    mpz_t y1;
    mpz_init_set_si(y1, y);
    mpz_add(result, x, y1);
    mpz_mod(result, result, PRIME);
    mpz_clear(y1);
}

void modMath::modAdd(mpz_t *result, mpz_t *x, long y, int size) {
    for (int i = 0; i < size; ++i)
        modAdd(result[i], x[i], y);
}

void modMath::modAdd(mpz_t *result, mpz_t *x, mpz_t y, int size) {
    for (int i = 0; i < size; i++)
        modAdd(result[i], x[i], y);
}

void modMath::modAdd(mpz_t *result, mpz_t *x, long *y, int size) {
    mpz_t *ytmp = (mpz_t *)malloc(sizeof(mpz_t) * size);
    for (int i = 0; i < size; i++)
        mpz_init_set_si(ytmp[i], y[i]);
    modAdd(result, x, ytmp, size);
    for (int i = 0; i < size; i++)
        mpz_clear(ytmp[i]);
}

void modMath::modAdd(mpz_t *result, mpz_t *x, int *y, int size) {
    mpz_t *ytmp = (mpz_t *)malloc(sizeof(mpz_t) * size);
    for (int i = 0; i < size; i++)
        mpz_init_set_si(ytmp[i], y[i]);
    modAdd(result, x, ytmp, size);
    for (int i = 0; i < size; i++)
        mpz_clear(ytmp[i]);
}

void modMath::modSub(mpz_t result, mpz_t x, mpz_t y) {
    mpz_sub(result, x, y);
    mpz_mod(result, result, PRIME);
}

void modMath::modSub(mpz_t *result, mpz_t *x, mpz_t *y, int size) {
    for (int i = 0; i < size; i++)
        modSub(result[i], x[i], y[i]);
}

void modMath::modSub(mpz_t result, mpz_t x, long y) {
    mpz_t y1;
    mpz_init_set_si(y1, y);
    mpz_sub(result, x, y1);
    mpz_mod(result, result, PRIME);
    mpz_clear(y1);
}

void modMath::modSub(mpz_t result, long x, mpz_t y) {
    mpz_t x1;
    mpz_init_set_si(x1, x);
    mpz_sub(result, x1, y);
    mpz_mod(result, result, PRIME);
    mpz_clear(x1);
}

void modMath::modSub(mpz_t *result, mpz_t *x, long y, int size) {
    for (int i = 0; i < size; ++i)
        modSub(result[i], x[i], y);
}

void modMath::modSub(mpz_t *result, long x, mpz_t *y, int size) {
    for (int i = 0; i < size; ++i)
        modSub(result[i], x, y[i]);
}

void modMath::modSub(mpz_t *result, mpz_t *x, mpz_t y, int size) {
    for (int i = 0; i < size; ++i)
        modSub(result[i], x[i], y);
}

void modMath::modSub(mpz_t *result, mpz_t x, mpz_t *y, int size) {
    for (int i = 0; i < size; ++i)
        modSub(result[i], x, y[i]);
}

void modMath::modSub(mpz_t *result, mpz_t *x, int *y, int size) {
    mpz_t *ytmp = (mpz_t *)malloc(sizeof(mpz_t) * size);
    for (int i = 0; i < size; i++)
        mpz_init_set_si(ytmp[i], y[i]);
    modSub(result, x, ytmp, size);
    for (int i = 0; i < size; i++)
        mpz_clear(ytmp[i]);
}
void modMath::modSub(mpz_t *result, int *x, mpz_t *y, int size) {
    mpz_t *xtmp = (mpz_t *)malloc(sizeof(mpz_t) * size);
    for (int i = 0; i < size; i++)
        mpz_init_set_si(xtmp[i], x[i]);
    modSub(result, xtmp, y, size);
    for (int i = 0; i < size; i++)
        mpz_clear(xtmp[i]);
}

void modMath::modPow2(mpz_t result, int exponent) {
    mpz_t value, base;
    mpz_init_set_ui(base, 2);
    mpz_init_set_si(value, exponent);
    // modAdd(value, value, (long)0); // assuming this just performs modular reduction, replaced with line below
    mpz_mod(value, value, PRIME);
    mpz_powm(result, base, value, PRIME);
    mpz_clear(value);
    mpz_clear(base);
}

void modMath::modPow2(mpz_t result, mpz_t exponent) {
    mpz_t value, base;
    mpz_init_set_ui(base, 2);
    mpz_init_set(value, exponent);
    // modAdd(value, value, (long)0); // assuming this just performs modular reduction, replaced with line below
    mpz_mod(value, value, PRIME);
    mpz_powm(result, base, value, PRIME);
    mpz_clear(value);
    mpz_clear(base);
}

void modMath::modPow2(mpz_t *result, int *exponent, int size) {
    // for (int i = 0; i < size; ++i)
    //     modPow(result[i], base[i], exponent);
    mpz_t value, base;
    mpz_init_set_ui(base, 2);

    for (int i = 0; i < size; ++i) {
        mpz_init_set_si(value, exponent[i]);
        mpz_mod(value, value, PRIME);
        mpz_powm(result[i], base, value, PRIME);
    }
    mpz_clear(value);
    mpz_clear(base);
}

void modMath::modPow2(mpz_t *result, mpz_t *exponent, int size) {
    // for (int i = 0; i < size; ++i)
    //     modPow(result[i], base[i], exponent);
    mpz_t value, base;
    mpz_init_set_ui(base, 2);

    for (int i = 0; i < size; ++i) {
        mpz_init_set(value, exponent[i]);
        mpz_mod(value, value, PRIME);
        mpz_powm(result[i], base, value, PRIME);
    }
    mpz_clear(value);
    mpz_clear(base);
}

void modMath::modPow(mpz_t result, mpz_t base, mpz_t exponent) {
    mpz_powm(result, base, exponent, PRIME);
}

void modMath::modPow(mpz_t *result, mpz_t *base, mpz_t *exponent, int size) {
    for (int i = 0; i < size; i++)
        mpz_powm(result[i], base[i], exponent[i], PRIME);
}

void modMath::modPow(mpz_t result, mpz_t base, long exponent) {
    mpz_t value;
    mpz_init_set_si(value, exponent);
    // modAdd(value, value, (long)0);
    mpz_mod(value, value, PRIME);
    mpz_powm(result, base, value, PRIME);
    mpz_clear(value);
}

void modMath::modPow(mpz_t *result, mpz_t *base, long exponent, int size) {
    // for (int i = 0; i < size; ++i)
    //     modPow(result[i], base[i], exponent);
    mpz_t value;
    mpz_init_set_si(value, exponent);
    mpz_mod(value, value, PRIME);
    for (int i = 0; i < size; ++i) {
        mpz_powm(result[i], base[i], value, PRIME);
    }
    mpz_clear(value);
}

void modMath::newModInv(mpz_t result, mpz_t value) {

    mpz_t temp1, temp2;
    mpz_init(temp1);
    mpz_init(temp2);
    // mpz_gcdext(result, temp1, temp2, value, PRIME);
    mpz_gcdext(temp1, result, temp2, value, PRIME);

    // gmp_printf("res %Zd \n", result);
    // gmp_printf("t1 %Zd \n", temp1);
    // gmp_printf("t2 %Zd \n", temp2);
    if (mpz_cmp_si(result, 0) > 0 or mpz_cmp_si(result, 0) == 0) {
        // std::cout<<"if\n";
        return;
    } else {

        // std::cout<<"else\n";
        mpz_add(result, result, PRIME);
        // mpz_mod(result, result, PRIME);
    }
}

void modMath::modInv(mpz_t result, mpz_t value) {
    mpz_t temp;
    mpz_init(temp);
    mpz_sub_ui(temp, PRIME, 2);
    modPow(result, value, temp);
    mpz_clear(temp);
}

void modMath::modInv(mpz_t *result, mpz_t *values, int size) {

    mpz_t temp;
    mpz_init(temp);
    mpz_sub_ui(temp, PRIME, 2);
    for (int i = 0; i < size; i++)
        modPow(result[i], values[i], temp);
    // modInv(result[i], values[i]); // highly inefficient
    mpz_clear(temp);
}

void modMath::modSqrt(mpz_t result, mpz_t x) {
    mpz_t temp;
    mpz_init(temp);
    mpz_add_ui(temp, PRIME, 1);
    mpz_div_ui(temp, temp, 4);
    modPow(result, x, temp);
    mpz_clear(temp);
}

void modMath::modSqrt(mpz_t *result, mpz_t *x, int size) {
    mpz_t *power = (mpz_t *)malloc(sizeof(mpz_t) * size);
    for (int i = 0; i < size; i++) {
        mpz_init(power[i]);
        mpz_add_ui(power[i], PRIME, 1);
        mpz_div_ui(power[i], power[i], 4);
    }
    modPow(result, x, power, size);
    for (int i = 0; i < size; i++)
        mpz_clear(power[i]);
}

void modMath::modSum(mpz_t result, mpz_t *x, int size) {
    mpz_set_ui(result, 0);
    for (int i = 0; i < size; ++i) {
        mpz_add(result, result, x[i]);
    }
    mpz_mod(result, result, PRIME);
}

void modMath::mod(mpz_t *result, mpz_t *a, mpz_t *m, int size) {
    mpz_t tmp;
    mpz_init(tmp);
    for (int i = 0; i < size; ++i) {
        mpz_init_set_ui(tmp, 0);
        mpz_add(tmp, a[i], m[i]);
        if (mpz_cmp(tmp, PRIME) > 0)
            mpz_sub(result[i], tmp, PRIME);
        else {
            mpz_mod(result[i], a[i], m[i]);
            mpz_mod(result[i], result[i], PRIME);
        }
    }
}

void modMath::mod(mpz_t *result, mpz_t *a, mpz_t m, int size) {
    mpz_t tmp;
    mpz_init(tmp);
    for (int i = 0; i < size; ++i) {
        mpz_init_set_ui(tmp, 0);
        mpz_add(tmp, a[i], m);
        if (mpz_cmp(tmp, PRIME) > 0)
            mpz_sub(result[i], tmp, PRIME);
        else {
            mpz_mod(result[i], a[i], m);
            mpz_mod(result[i], result[i], PRIME);
        }
    }
}

void modMath::copy(mpz_t *src, mpz_t *des, int size) {
    for (int i = 0; i < size; i++)
        mpz_set(des[i], src[i]);
}
