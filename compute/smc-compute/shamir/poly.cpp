
#include "poly.hpp"
#include <algorithm>
#include <cstdio>
#include <functional>
#include <iostream>

int poly::computeTrueLen(vector<int> A, int A_len) {
    for (int i = A_len - 1; i >= 0; i--) {
        if (A[i] != 0) {
            return i + 1;
        }
    }
    return 0;
}

int poly::computeTrueLen(mpz_t *A, int A_len) {
    for (int i = A_len - 1; i >= 0; i--) {
        if (mpz_cmp_si(A[i], 0) != 0) {
            return i + 1;
        }
    }
    return 0;
}

// checks if there are any trailing zero coeffs
// shortens coeffs if necessary, updates the coef_sz and degree
void poly::truncateZeros() {

    int new_coef_sz = computeTrueLen(this->coeffs, this->coef_sz);

    if (new_coef_sz < this->coef_sz) {

        mpz_t *new_coeffs = (mpz_t *)malloc(sizeof(mpz_t) * new_coef_sz);

        for (int i = 0; i < new_coef_sz; i++) {
            mpz_init(new_coeffs[i]);
        }
        copy(this->coeffs, new_coeffs, new_coef_sz);
        freeCoeffs();
        this->coeffs = new_coeffs;
        this->coef_sz = new_coef_sz;
        this->degree = new_coef_sz - 1;
    }
}

poly::poly(mpz_t _PRIME, mpz_t *_coef, int _coef_sz) : modMath(_PRIME) {

    // printf("CONSTRUCTOR CALLED\n");
    coef_sz = computeTrueLen(_coef, _coef_sz);
    // std::cout << "coef_sz = " << coef_sz << std::endl;
    if (coef_sz > 0) {

        coeffs = (mpz_t *)malloc(sizeof(mpz_t) * coef_sz);
        for (int i = 0; i < coef_sz; i++) {
            mpz_init(coeffs[i]);
        }
        copy(_coef, coeffs, coef_sz); // need to make a copy here so this has ownership

        // coeffs = _coef;
        degree = coef_sz - 1;
    } else {
        degree = -1;
    }
}

// init zero poly
poly::poly(mpz_t _PRIME) : modMath(_PRIME) {
    coef_sz = 0;
    degree = -1;
}

poly::poly(mpz_t _PRIME, vector<int> _coef) : modMath(_PRIME) {

    // printf("CONSTRUCTOR CALLED\n");
    coef_sz = computeTrueLen(_coef, _coef.size());

    // std::cout << "v coef_sz = " << coef_sz << std::endl;
    if (coef_sz > 0) {

        coeffs = (mpz_t *)malloc(sizeof(mpz_t) * coef_sz);
        for (size_t i = 0; i < coef_sz; i++) {
            mpz_init_set_si(coeffs[i], _coef.at(i));
        }
        degree = coef_sz - 1;
    } else {
        degree = -1;
    }
    // mod(coeffs, coeffs, coef_sz);
}

void poly::print(std::string name) {

    cout << name << ": ";
    for (size_t j = 0; j < this->coef_sz; j++) {
        gmp_printf("%Zd, ", this->coeffs[j]);
    }
    printf("\n");
}

void poly::updatePoly(mpz_t *new_coeff, int new_coef_sz) {
    freeCoeffs();

    coef_sz = computeTrueLen(new_coeff, new_coef_sz);

    if (coef_sz > 0) {
        coeffs = (mpz_t *)malloc(sizeof(mpz_t) * coef_sz);

        for (int i = 0; i < coef_sz; i++) {
            mpz_init(coeffs[i]);
        }
        copy(new_coeff, coeffs, coef_sz); // need to make a copy here so this has ownership
        degree = coef_sz - 1;
    } else {
        degree = -1;
    }
}

poly::~poly() {
    // printf("DESTRUCTOR CALLED\n");
    freeCoeffs();
}

// calls mpz_clear on coeffs and deletes the memory
void poly::freeCoeffs() {
    if (degree != -1) {
        for (int i = 0; i < coef_sz; i++) {
            mpz_clear(coeffs[i]);
        }

    free(coeffs);
    }
}

const std::size_t poly::deg() const {
    return degree;
}

void poly::plus(poly &other) {
    if (this->deg() == -1 && other.deg() != -1) {

        // need to make a copy here so this has ownership
        // we are effectively creating a copy

        // coeffs was never initalized
        this->coef_sz = other.coef_sz;

        coeffs = (mpz_t *)malloc(sizeof(mpz_t) * other.coef_sz);

        for (int i = 0; i < coef_sz; i++) {
            mpz_init(coeffs[i]);
        }
        copy(other.coeffs, this->coeffs, this->coef_sz);
        this->degree = this->coef_sz - 1;
        return;
    } else if (this->deg() != -1 && other.deg() == -1) {
        return;
        // do nothing, since the other poly is the zero poly
    }
    if (other.deg() > this->deg()) {
        // need to reallocate our coeffs to be larger

        // mpz_t *new_coeffs = new mpz_t[other.coef_sz];
        mpz_t *new_coeffs = (mpz_t *)malloc(sizeof(mpz_t) * other.coef_sz);

        for (int i = 0; i < other.coef_sz; i++) {
            mpz_init(new_coeffs[i]);
        }
        copy(this->coeffs, new_coeffs, this->coef_sz);
        freeCoeffs();
        this->coeffs = new_coeffs;
        this->coef_sz = other.coef_sz;
        this->degree = other.degree;
    }

    modAdd(this->coeffs, this->coeffs, other.coeffs, this->coef_sz);

    // for (size_t i = 0; i < this->coef_sz; i++) {
    // gmp_printf("plus this (%i, %Zd) \n", i, this->coeffs[i]);
    // }
    truncateZeros();
}

void poly::sub(poly &other) {
    if (this->deg() == -1 && other.deg() != -1) {

        // need to make a copy here so this has ownership
        // we are effectively creating a copy

        // coeffs was never initalized
        this->coef_sz = other.coef_sz;

        this->coeffs = (mpz_t *)malloc(sizeof(mpz_t) * other.coef_sz);
        // this->coeffs = new mpz_t[other.coef_sz];

        for (int i = 0; i < other.coef_sz; i++) {
            mpz_init(coeffs[i]);
            // }
            // for (int i = 0; i < other.coef_sz; i++) {
            modSub(this->coeffs[i], long(0), other.coeffs[i]);
        }

        // copy(other.coeffs, this->coeffs, this->coef_sz);
        this->degree = this->coef_sz - 1;
        return;
    } else if (this->deg() != -1 && other.deg() == -1) {
        return;
        // do nothing, since the other poly is the zero poly
    }
    if (other.deg() > this->deg()) {
        // need to reallocate our coeffs to be larger


        mpz_t *new_coeffs = (mpz_t *)malloc(sizeof(mpz_t) * other.coef_sz);
        // mpz_t *new_coeffs = new mpz_t[other.coef_sz];

        for (int i = 0; i < other.coef_sz; i++) {
            mpz_init(new_coeffs[i]);
        }
        copy(this->coeffs, new_coeffs, this->coef_sz);
        freeCoeffs();
        this->coeffs = new_coeffs;
        this->coef_sz = other.coef_sz;
        this->degree = other.degree;
    }

    modSub(this->coeffs, this->coeffs, other.coeffs, this->coef_sz);
    truncateZeros();
}

void poly::multScalar(mpz_t val) {
    if (this->deg() == -1) {
        // either of the polys are the zero polynomial
        // cout << "at least one operand was the zero poly" << endl;
        if (coef_sz > 0) {
            freeCoeffs();
        }
        this->coef_sz = 0;
        this->degree = -1;
        return;
    }
    // need to reallocate our coeffs to be larger

    // size_t new_coeff_sz = this->degree + other.degree + 1;
    // mpz_t *new_coeffs = new mpz_t[new_coeff_sz];
    // for (int i = 0; i < new_coeff_sz; i++) {
    //     mpz_init(new_coeffs[i]);
    // }

    // mpz_t temp;
    // mpz_init(temp);
    for (int i = 0; i < this->coef_sz; i++) {
        // gmp_printf("coef[%i] %Zd, val %Zd \n", i, this->coeffs[i], val);
        // gmp_printf("res[%i] %Zd, \n", i, this->coeffs[i]);
    }
    modMul(this->coeffs, this->coeffs, val, this->coef_sz);
    truncateZeros();
}

void poly::mult(poly &other) {
    if (this->deg() == -1 || other.deg() == -1) {
        // either of the polys are the zero polynomial

        // cout << "at least one operand was the zero poly" << endl;
        if (coef_sz > 0) {
            freeCoeffs();
        }
        this->coef_sz = 0;
        this->degree = -1;
        return;
    }
    // need to reallocate our coeffs to be larger

    size_t new_coeff_sz = this->degree + other.degree + 1;
    mpz_t *new_coeffs = new mpz_t[new_coeff_sz];
    for (int i = 0; i < new_coeff_sz; i++) {
        mpz_init(new_coeffs[i]);
    }

    mpz_t temp;
    mpz_init(temp);
    for (int i = 0; i < this->coef_sz; i++) {
        for (int j = 0; j < other.coef_sz; j++) {
            modMul(temp, this->coeffs[i], other.coeffs[j]);
            modAdd(new_coeffs[i + j], new_coeffs[i + j], temp);
        }
    }
    freeCoeffs();

    // need to re-assign new_coeffs to this->coeffs
    this->coeffs = new_coeffs;
    this->coef_sz = new_coeff_sz;
    this->degree = new_coeff_sz - 1;

    truncateZeros();
}

// computes this (A) /  other (B), stores quotient poly in quotient and remainder in remainder
void poly::div_mod(poly &quotient, poly &remainder, poly &other) {
    mpz_t temp;
    mpz_t inv;
    mpz_init(temp);
    mpz_init(inv);

    // std::cout << "coef sz" << this->coef_sz << "\n";
    mpz_t *Q = (mpz_t *)malloc(sizeof(mpz_t) * this->coef_sz);
    mpz_t *R = (mpz_t *)malloc(sizeof(mpz_t) * this->coef_sz);
    // mpz_t *Q = new mpz_t[this->coef_sz];
    // mpz_t *R = new mpz_t[this->coef_sz];
    for (int i = 0; i < this->coef_sz; i++) {
        mpz_init(Q[i]);                      // init to zero
        mpz_init_set(R[i], this->coeffs[i]); // init to numerator
    }

    // ensure there is no trailing zeros in other so we can get the coeff of highest order term
    // but we shouldn't ever need to? when we initalize we make sure its not the case
    // other.truncateZeros();

    // gmp_printf("inv %Zd \n", other.coeffs[other.coef_sz - 1]);
    modInv(inv, other.coeffs[other.coef_sz - 1]);
    // gmp_printf("inv %Zd \n", inv);

    int len = this->coef_sz - other.coef_sz + 1;

    for (int i = len - 1; i >= 0; i--) {
        modMul(Q[i], R[i + other.coef_sz - 1], inv);
        for (int j = 0; j < other.coef_sz; j++) {
            modMul(temp, Q[i], other.coeffs[j]);
            modSub(R[i + j], R[i + j], temp);
        }
    }
    quotient.updatePoly(Q, this->coef_sz);
    remainder.updatePoly(R, this->coef_sz);

    // Q_size = computePolyLen(quotient, A_size);
    // R_size = computePolyLen(remainder, A_size);
    // compute output length of Q and R (stripping trailing zeros)

    mpz_clear(temp);
    mpz_clear(inv);
}

void generateLagrangeCoeff(std::vector<poly *> &lagrange_poly, vector<int> &points, mpz_t modulus) {

    modMath math = modMath(modulus);
    mpz_t nom, denom, t1, t2, temp;
    mpz_init(nom);
    mpz_init(denom);
    mpz_init(t1);
    mpz_init(t2);
    mpz_init(temp);
    int xi, xj;
    poly temp_poly = poly(modulus, {1, 1});
    for (size_t i = 0; i < points.size(); i++) {

        mpz_set_ui(denom, 1);

        mpz_set_ui(t2, points.at(i)); // xi
        poly *numerator = new poly(modulus, {1});
        for (size_t j = 0; j < points.size(); j++) {
            if (i != j) {
                // xj = points.at(j);
                mpz_set_ui(t1, points.at(j));                  // xj
                math.modSub(temp_poly.coeffs[0], long(0), t1); //  -xj
                numerator->mult(temp_poly);

                math.modSub(temp, t2, t1); //  xi - xj
                math.modMul(denom, denom, temp);
            }
        }

        // for (size_t i = 0; i < numerator->coef_sz; i++) {
        //     gmp_printf("numerator (%i, %Zd) \n", i, numerator->coeffs[i]);
        // }
        //
        // gmp_printf("denom %Zd \n", denom);
        math.modInv(temp, denom);
        math.modMul(numerator->coeffs, numerator->coeffs, temp, numerator->coef_sz);

        // for (size_t i = 0; i < numerator->coef_sz; i++) {
        //     gmp_printf("num after (%i, %Zd) \n", i, numerator->coeffs[i]);
        // }
        // printf("\n");
        // printf("\n");
        lagrange_poly.push_back(numerator);
        // std::cout << "hi\n";
        // for (int l = 0; l < peers; l++) {
        //     modMul(lagrangePolysAll[i][l], num_poly[l], temp);
        // }
    }

    // delete temp_poly;
    // std::cout << "end\n";
}

void interpolate(poly &result, vector<int> &points, mpz_t *shares, int numShares, mpz_t modulus) {
    modMath math = modMath(modulus);
    std::vector<poly *> lagrange_poly = {};
    generateLagrangeCoeff(lagrange_poly, points, modulus);
    if (points.size() != numShares) {
        cout << "points and numShares must be equal\n";
        return;
    }

    // for (auto row : lagrange_poly) {
    //     printf("ls ");
    //     for (size_t i = 0; i < row->coef_sz; i++) {
    //         gmp_printf("%Zd, ", row->coeffs[i]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    for (size_t i = 0; i < numShares; i++) {
        // math.modMul(lagrange_poly.at(i)->coeffs, lagrange_poly.at(i)->coeffs, shares[i], lagrange_poly.at(i)->coef_sz);
        lagrange_poly.at(i)->multScalar(shares[i]);
        // gmp_printf("share (%i, %Zd) \n", i, shares[i]);

        // for (size_t j = 0; j < lagrange_poly.at(i)->coef_sz; j++) {
        //     gmp_printf("term (%j, %Zd) \n", j,  lagrange_poly.at(i)->coeffs[j]);
        // }
        // printf("\n");

        // printf("term ");
        // for (size_t j = 0; j < lagrange_poly.at(i)->coef_sz; j++) {
        //     gmp_printf("%Zd, ",  lagrange_poly.at(i)->coeffs[j]);
        // }
        // printf("\n");

        result.plus(*lagrange_poly.at(i));

        // printf("res ");
        // for (size_t j = 0; j < result.coef_sz; j++) {
        //     gmp_printf("%Zd, ",  result.coeffs[j]);
        // }
        // printf("\n");
        // printf("\n");

        // std::cout <<"updated "<<result.coef_sz << endl;
        // for (size_t i = 0; i < result.coef_sz; i++) {
        //     gmp_printf("res (%i, %Zd) \n", i, result.coeffs[i]);
        // }

        // printf("\n");
    }

    // cleanup
    for (auto p : lagrange_poly)
        delete p;
    lagrange_poly.clear();
}

void RS_decode(poly &result, poly &error_loc, vector<int> &points, mpz_t *shares, int numShares, int max_degree, int max_error_count, mpz_t modulus) {

    // std::cout << "hi\n";
    poly H = poly(modulus);
    interpolate(H, points, shares, numShares, modulus);

    modMath math = modMath(modulus);
    mpz_t t1, t2;
    mpz_init(t1);
    mpz_init(t2);

    poly temp_poly = poly(modulus, {1, 1});
    poly F = poly(modulus, {1});

    for (size_t i = 0; i < points.size(); i++) {

        mpz_set_ui(t1, points.at(i));                  // xj
        math.modSub(temp_poly.coeffs[0], long(0), t1); //  -xj
        F.mult(temp_poly);
    }

    poly R0 = F;
    poly R1 = H;
    poly R2 = poly(modulus);

    poly S0 = poly(modulus, {1});
    poly S1 = poly(modulus, {});
    poly T0 = poly(modulus, {});
    poly T1 = poly(modulus, {1});

    poly R = poly(modulus);
    poly Q = poly(modulus);

    poly temp1 = poly(modulus);
    poly temp2 = poly(modulus);

    poly G = poly(modulus);
    poly remainder = poly(modulus);
    while (true) {
        R0.div_mod(Q, R2, R1);

        if (R0.degree < (max_error_count + max_degree)) {
            R0.div_mod(G, remainder, T0);

            // remainder is the zero polynomial
            if (remainder.degree == -1) {
                result = G;
                error_loc = T0;

                return;
            } else {
                cout << "too many errors encountered, aborting..." << endl;
                return;
            }
        }

        R0 = R1;
        R1 = R2;

        temp1 = S1;
        temp1.mult(Q);
        temp2 = S0;
        temp2.sub(temp1);

        S0 = S1;
        S1 = temp2;

        temp1 = T1;

        temp2 = T0;
        temp2.sub(temp1);
        T0 = T1;
        T1 = temp2;
    }
}
