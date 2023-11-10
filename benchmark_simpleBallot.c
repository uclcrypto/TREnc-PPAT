#include <stdio.h> 
#include <stdlib.h> 

#include "source/TREnc_PPAT_simple_ballots.c"

int main()
{
    int i, iterations;
    clock_t start;
    double elapsed;
    char pr[10];
    unsigned long ran;

    csprng RNG;
    BIG_BN462 m, decrypted_m;
    sk_simple_ballot SK;
    pk_simple_ballot PK;
    simple_ballot ballot;

    BIG_464_60_one(m);

    printf("\nTesting/Timing the scheme for simple ballots\n");

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
        keygen(&RNG, &SK, &PK);
    }
    elapsed = (clock() - start) / (double)CLOCKS_PER_SEC;
    elapsed = 1000.0 * elapsed/iterations;
    printf("KEYGEN RUNNING TIME             - %8d iterations  ", iterations);
    printf(" %8.2lf ms per iteration\n", elapsed);

    start = clock();
    for (i = 0; i < iterations; i++){
        encrypt(&RNG, m, PK, &ballot);
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
   
    decrypt(decrypted_m, SK, PK, ballot);
    if (BIG_464_60_comp(decrypted_m, m) == 0)
        printf("DECRYPTION IS CORRECT!\n");   
    else
        printf("DECRYPTION IS NOT CORRECT!\n");

  
/* 
    // COMPUTE FOR K BITS
    int k =10;
    int randomBit = rand() % 2;
    BIG_BN462 array_m[k];
    for (i = 0; i < k; i++){
        if (randomBit ==0){
            BIG_464_60_zero(array_m[i]);
        }
        else {
            BIG_464_60_one(array_m[i]);
        }
    }

    start = clock();
    for (i = 0; i < k; i++){
        enc(&RNG, array_m[i]);
    }
    elapsed = (clock() - start) / (double)CLOCKS_PER_SEC;
    elapsed = 1000.0 * elapsed;
    printf("ENCRYPTION RUNNING TIME FOR K %8.2lf ms \n", elapsed);
    */
   return 0;
}
