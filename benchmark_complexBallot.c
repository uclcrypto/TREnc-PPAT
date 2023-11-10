#include <stdio.h> 
#include <stdlib.h> 

#include "source/TREnc_PPAT_complex_ballots.c"

int main()
{
    int i, iterations;
    clock_t start;
    double elapsed;
    char pr[10];
    unsigned long ran;

    bool is_valid;
    csprng RNG;
    BIG_BN462  m, decrypted_m, order_group;
    ECP_BN462 M, decrypted_M;
    sk_complex_ballot SK;
    pk_complex_ballot PK;
    complex_ballot ballot;

    // Pick a random message m;
    // and convert it to a point in BN462 curve.
    BIG_BN462_rcopy(order_group, CURVE_Order_BN462);
    do{
    BIG_BN462_randtrunc(m, order_group, 2 * CURVE_SECURITY_BN462, &RNG);
    is_valid = ECP_BN462_setx(&M, m, 1);
    } while (is_valid == false);
    
    printf("\nTesting/Timing the scheme for complex ballots\n");

    time((time_t *)&ran);
    pr[0] = ran; 
    pr[1] = ran >> 8;
    pr[2] = ran >> 16;
    pr[3] = ran >> 24;
    for (i = 4; i < 10; i++){
        pr[i] = i;
    } 
    RAND_seed(&RNG, 10, pr);

    iterations = 10;
    start = clock();
    for (i = 0; i < iterations; i++){
        generate_keys(&RNG, &SK, &PK);
    }
    elapsed = (clock() - start) / (double)CLOCKS_PER_SEC;
    elapsed = 1000.0 * elapsed/iterations;
    printf("KEYGEN RUNNING TIME             - %8d iterations  ", iterations);
    printf(" %8.2lf ms per iteration\n", elapsed);

    start = clock();
    for (i = 0; i < iterations; i++){
        encrypt(&RNG, M, PK, &ballot);
    }
    
    elapsed = (clock() - start) / (double)CLOCKS_PER_SEC;
    elapsed = 1000.0 * elapsed/iterations;
    printf("ENCRYPTION RUNNING TIME         - %8d iterations  ", iterations);
    printf(" %8.2lf ms per iteration\n", elapsed);

    start = clock();
    for (i = 0; i < iterations; i++){
        verify(PK, ballot);
    }
    elapsed = (clock() - start) / (double)CLOCKS_PER_SEC;
    elapsed = 1000.0 * elapsed/iterations;
    printf("VERIFICATION RUNNING TIME       - %8d iterations  ", iterations);
    printf(" %8.2lf ms per iteration\n", elapsed);

    decrypt(&decrypted_M, SK, PK, ballot);
    ECP_BN462_get(decrypted_m, decrypted_m, &decrypted_M);    
    if (BIG_464_60_comp(decrypted_m, m) == 0)
        printf("DECRYPTION IS CORRECT!\n");
    else
        printf("DECRYPTION IS NOT CORRECT!\n");

    return 0;
}
