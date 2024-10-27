#include <float.h>
#include <iomanip>
#include <limits.h>

extern "C" int ort_initialize(int *, char ***);
extern "C" void ort_finalize(int);

#include "smc-compute/SMC_Utils.h"
#include <gmp.h>

SMC_Utils *__s;

std::vector<int> seed_map;
// bool ERROR_FLAG = false;

int __original_main(int _argc_ignored, char **_argv_ignored) {

    int L = 4096;
    priv_int bprime;

    ss_init(bprime);
    priv_int *aprime;
    priv_int *sprime;
    sprime = (priv_int *)malloc(sizeof(priv_int) * (L));
    aprime = (priv_int *)malloc(sizeof(priv_int) * (L));
    for (int _picco_i = 0; _picco_i < L; _picco_i++) {
        ss_init(sprime[_picco_i]);
        ss_init(aprime[_picco_i]);
    }

    __s->smc_input(1, sprime, L, "int", -1);
    __s->smc_input(2, aprime, L, "int", -1);
    __s->smc_input(2, &bprime, "int", -1);

    for (int ERROR_FLAG = 0; ERROR_FLAG < 2; ERROR_FLAG++) {
        int iterations = 1;
        double offline_time = 0.0, online_time = 0.0;
        std::cout << "----------------------------------------" << std::endl;

        for (size_t i = 0; i < iterations; i++) {
            __s->thresholdDecryption(sprime, aprime, bprime, L, -1, offline_time, online_time, bool(ERROR_FLAG));
        }

        std::cout << std::boolalpha;
        std::cout << "Errors introduced? " << bool(ERROR_FLAG) << std::endl;
        std::cout << "Offline (PRSS):         " << offline_time / iterations << " ms " << std::endl;
        std::cout << "Online (v, robustOpen): " << online_time / iterations << " ms " << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "Total runtime:          " << (offline_time + online_time) / iterations << " ms " << std::endl;
    }
    for (int _picco_i = 0; _picco_i < L; _picco_i++) {
        ss_clear(sprime[_picco_i]);
        ss_clear(aprime[_picco_i]);
    }
    free(sprime);
    free(aprime);
    ss_clear(bprime);

    return 0;
}

/* smc-compiler generated main() */
int main(int argc, char **argv) {

    if (argc < 6) {
        fprintf(stderr, "Incorrect input parameters\n");
        fprintf(stderr, "Usage: <id> <runtime-config> <privatekey-filename> <numParties> <sprime_share> <ciphertext>\n");
        exit(1);
    }

    int numParties = atoi(argv[4]);
    int threshold;

    std::string IO_files[2];
    for (int i = 0; i < 2; i++) {
        IO_files[i] = argv[5 + i];
        std::cout << IO_files[i] << std::endl;
    }

    switch (numParties) {
    case 5:
        seed_map = {7, 11};
        threshold = 2;
        break;
    case 11:
        seed_map = {63, 95, 111, 119, 123, 125, 159, 175, 183, 187, 189, 207, 215, 219, 221, 231, 235, 237, 243, 245, 249, 303, 311, 315, 317, 335, 343, 347, 349, 359, 363, 365, 371, 373, 411, 413, 423, 427, 429, 437, 469, 683};
        threshold = 5;
        break;
    default:
        std::cerr << "wrong number of parties!" << std::endl;
        exit(1);
        break;
    }

        std::cout << "numParties "<< numParties << std::endl;
        std::cout << "threshold "<< threshold << std::endl;
    __s = new SMC_Utils(atoi(argv[1]), argv[2], argv[3], 2, 0, IO_files, numParties, threshold, 128, "200393528695012829568562844035540003081", seed_map, 1);
    // __s = new SMC_Utils(atoi(argv[1]), argv[2], argv[3], 2, 0, IO_files, numParties, threshold, 128, "170141183460469231731687303715884105851", seed_map, 1);

    struct timeval tv1;
    struct timeval tv2;
    int _xval = 0;

    gettimeofday(&tv1, NULL);
    _xval = (int)__original_main(argc, argv);
    gettimeofday(&tv2, NULL);
    return (_xval);
}
