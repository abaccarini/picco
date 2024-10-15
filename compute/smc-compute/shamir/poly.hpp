#pragma once

#include "modMath.hpp"
#include "SecretShare.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

class poly : public modMath {
public:
    poly(mpz_t, mpz_t *, int);
    poly(mpz_t, vector<int>);
    poly(mpz_t);
    int computeTrueLen(mpz_t *, int);
    int computeTrueLen(vector<int>, int);
    ~poly();

    poly(poly &t) : modMath(t) {
        this->coef_sz = t.coef_sz;
        coeffs = new mpz_t[this->coef_sz];
        for (int i = 0; i < this->coef_sz; i++) {
            mpz_init(coeffs[i]);
        }
        copy(t.coeffs, this->coeffs, this->coef_sz); // need to make a copy here so this has ownership

        this->degree = coef_sz - 1;
    };

    poly &operator=(const poly &other) {

        this->coef_sz = other.coef_sz;
        coeffs = new mpz_t[this->coef_sz];
        for (int i = 0; i < this->coef_sz; i++) {
            mpz_init(coeffs[i]);
        }
        copy(other.coeffs, this->coeffs, this->coef_sz); // need to make a copy here so this has ownership

        this->degree = coef_sz - 1;
        return *this;
    }
    void updatePoly(mpz_t *new_coeff, int new_coef_sz);

    void freeCoeffs();
    const std::size_t deg() const;

    void print(std::string name);
    void plus(poly &);
    void sub(poly &);
    void mult(poly &);

    void multScalar(mpz_t val);
    void div_mod(poly &quotient, poly &remainder, poly &other);

    void truncateZeros();

    mpz_t *coeffs;
    int coef_sz;
    int degree; // the order
    // bool zero_poly;
private:
};

void generateLagrangeCoeff(std::vector<poly *> &lagrange_poly, vector<int> &points, mpz_t modulus);

void interpolate(poly &result, vector<int> &points, mpz_t *shares, int numShares, mpz_t modulus);

void RS_decode(poly &result, poly &error_loc, vector<int> &points, mpz_t *shares, int numShares, int max_degree, int max_error_count, mpz_t modulus);
