#include <float.h>
#include <limits.h>

extern "C" int ort_initialize(int *, char ***);
extern "C" void ort_finalize(int);

#include "smc-compute/SMC_Utils.h"
#include <gmp.h>

SMC_Utils *__s;

std::vector<int> seed_map = {7, 11};

int __original_main(int _argc_ignored, char **_argv_ignored) {

    int L = 3;

    priv_int bprime;
    priv_int *aprime;
    priv_int *sprime;
    sprime = (priv_int *)malloc(sizeof(priv_int) * (L));
    aprime = (priv_int *)malloc(sizeof(priv_int) * (L));
    for (int _picco_i = 0; _picco_i < L; _picco_i++) {
        ss_init(sprime[_picco_i]);
        ss_init(aprime[_picco_i]);
    }

    __s->smc_input(1, aprime, L, "int", -1);
    __s->smc_input(1, sprime, L, "int", -1);
    __s->smc_input(1, &bprime, "int", -1);

    __s->thresholdDecryption(sprime, aprime, bprime, L, -1);

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

    __s = new SMC_Utils(atoi(argv[1]), argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), IO_files, 5, 2, 128, "170141183460469231731687303715884105851", seed_map, 1);

    struct timeval tv1;
    struct timeval tv2;
    int _xval = 0;

    gettimeofday(&tv1, NULL);

    _xval = (int)__original_main(argc, argv);
    gettimeofday(&tv2, NULL);
    std::cout << "Time: " << __s->time_diff(&tv1, &tv2) << " seconds " << std::endl;
    return (_xval);
}
