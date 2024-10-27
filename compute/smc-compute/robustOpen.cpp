
#include "robustOpen.h"

// Robust opening based onReed-Solomon decoding
void RobustOpen(mpz_t result, mpz_t var, bool error_flag, int threadID, NodeNetwork nodeNet, SecretShare *ss) {
    uint peers = ss->getPeers();

    mpz_t modulus;
    mpz_init(modulus);
    ss->getFieldSize(modulus);

    mpz_t *shares = (mpz_t *)malloc(sizeof(mpz_t) * (peers));

    for (size_t i = 0; i < peers; i++) {
        mpz_init(shares[i]);
    }

    mpz_t *data = (mpz_t *)malloc(sizeof(mpz_t) * 1);
    // mpz_t *results = (mpz_t *)malloc(sizeof(mpz_t) * 1);
    mpz_init(data[0]);
    // mpz_init(results[0]);
    mpz_set(data[0], var);

    mpz_t **buffer = (mpz_t **)malloc(sizeof(mpz_t *) * (peers));
    for (int i = 0; i < (peers); i++) {
        buffer[i] = (mpz_t *)malloc(sizeof(mpz_t) * 1);
        for (int j = 0; j < 1; j++)
            mpz_init(buffer[i][j]);
    }
    nodeNet.broadcastToPeers(data, 1, buffer, threadID);
    for (size_t i = 0; i < peers; i++) {
        mpz_set(shares[i], buffer[i][0]); // just moving buffer into shares for consistency
        // gmp_printf("shares[%i] %Zd \n", i, shares[i]);
    }

    std::vector<int> points(peers);
    std::iota(std::begin(points), std::end(points), 1);
    // for (auto i : points) {
    //     std::cout << "point i  " << i << endl;
    // }
    int MAX_MANIPULATED = ss->getThreshold() - 1;
    int MAX_DEGREE = ss->getThreshold() + 1;

    if (error_flag) {
        // std::cout << "MAX_manipulated "<<MAX_MANIPULATED<<endl;
        for (size_t i = 0; i < MAX_MANIPULATED; i++) {
            // randomizing the maximum allowable shares to simulate a malicious party
            mpz_set_ui(shares[i], 0);
        }
    }

    poly out_poly = poly(modulus);
    poly error_loc = poly(modulus);


    // for (size_t i = 0; i < peers; i++) {
    //     gmp_printf("shares[%i] %Zd \n", i, shares[i]);
    // }
    RS_decode(out_poly, error_loc, points, shares, points.size(), MAX_DEGREE, MAX_MANIPULATED, modulus);

    // gmp_printf("result %Zd \n", out_poly.coeffs[0]);
    // mpz_set(result, out_poly.coeffs[0]); // segfault?

    mpz_clear(data[0]);
    free(data);

    // freeing
    for (int i = 0; i < (peers); i++) {
        for (int j = 0; j < 1; j++)
            mpz_clear(buffer[i][j]);
        free(buffer[i]);
    }
    free(buffer);
    // delete out_poly;
    // delete error_loc;
}
