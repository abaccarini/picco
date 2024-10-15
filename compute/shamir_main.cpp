#include <float.h>
#include <limits.h>
#include <iomanip>

extern "C" int ort_initialize(int *, char ***);
extern "C" void ort_finalize(int);

#include "smc-compute/SMC_Utils.h"
#include <gmp.h>

SMC_Utils *__s;

// std::vector<int> seed_map = {7, 11};
std::vector<int> seed_map = {63, 95, 111, 119, 123, 125, 159, 175, 183, 187, 189, 207, 215, 219, 221, 231, 235, 237, 243, 245, 249, 303, 311, 315, 317, 335, 343, 347, 349, 359, 363, 365, 371, 373, 411, 413, 423, 427, 429, 437, 469, 683};

bool ERROR_FLAG = true;

int __original_main(int _argc_ignored, char **_argv_ignored) {

    int L = 1024;

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

    int iterations = 1000;
    double offline_time = 0.0, online_time = 0.0;
    std::cout << "----------------------------------------" << std::endl;

    for (size_t i = 0; i < iterations; i++) {
        __s->thresholdDecryption(sprime, aprime, bprime, L, -1, offline_time, online_time, ERROR_FLAG);
    }

    std::cout << std::boolalpha;
    std::cout << "Errors introduced? " << ERROR_FLAG << std::endl;
    std::cout << "Offline (PRSS):         " << offline_time / iterations << " ms " << std::endl;
    std::cout << "Online (v, robustOpen): " << online_time / iterations << " ms " << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Total runtime:          " << (offline_time + online_time) / iterations << " ms " << std::endl;

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

    if (argc < 8) {
        fprintf(stderr, "Incorrect input parameters\n");
        fprintf(stderr, "Usage: <id> <runtime-config> <privatekey-filename> <number-of-input-parties> <number-of-output-parties> <input-share> <output>\n");
        exit(1);
    }

    std::string IO_files[atoi(argv[4]) + atoi(argv[5])];
    for (int i = 0; i < argc - 6; i++)
        IO_files[i] = argv[6 + i];

    // __s = new SMC_Utils(atoi(argv[1]), argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), IO_files, 5, 2, 128, "170141183460469231731687303715884105851", seed_map, 1);
    __s = new SMC_Utils(atoi(argv[1]), argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), IO_files, 11, 5, 128, "170141183460469231731687303715884105851", seed_map, 1);

    struct timeval tv1;
    struct timeval tv2;
    int _xval = 0;

    gettimeofday(&tv1, NULL);

    _xval = (int)__original_main(argc, argv);
    gettimeofday(&tv2, NULL);
    // std::cout << "Time: " << __s->time_diff(&tv1, &tv2) << " seconds " << std::endl;
    return (_xval);
}
