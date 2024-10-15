
#include "robustOpen.h"




// Robust opening based onReed-Solomon decoding
void RobustOpen(mpz_t result, mpz_t var, int threadID, NodeNetwork nodeNet, SecretShare *ss) {
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
    }
    std::vector<int> points(peers);                     
    std::iota(std::begin(points), std::end(points), 1); 



    int MAX_MANIPULATED = ss->getThreshold()-1;
    int MAX_DEGREE = ss->getThreshold()+1;
    poly *out_poly = new poly(modulus);
    poly *error_loc = new poly(modulus);
    RS_decode(*out_poly, *error_loc, points, shares, points.size(), MAX_DEGREE , MAX_MANIPULATED, modulus);

    mpz_set(result, out_poly->coeffs[0]);

    mpz_clear(data[0]);
    free(data);

    // freeing
    for (int i = 0; i < (peers); i++) {
        for (int j = 0; j < 1; j++)
            mpz_clear(buffer[i][j]);
        free(buffer[i]);
    }
    free(buffer);
    delete out_poly;
    delete error_loc;

}

