#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "tools.h"

void generate_keys(csprng *RNG, sk_complex_ballot *SK, pk_complex_ballot *PK){
    int i = 0, j = 0, len_sig = 2, len_sk = 3;
    ECP_BN462 gen1, g, h, g1, h1, temp1, f, crs_u[2][2], crs_u_prime[2][2], pk_sxdh[len_sig], temp_i[len_sig];
    ECP2_BN462  crs_v[2][2], gen2, g_hat, h_hat, f1_hat, f2_hat, temp2;
    BIG_BN462 r, x1, x2, r_u[2][2], r_u_prime[2][2], r_v[2][2], a[len_sk], b[len_sk], chi_i[len_sig], gamma_i[len_sig];

    ECP_BN462_generator(&gen1);
    ECP2_BN462_generator(&gen2);
    BIG_BN462_rcopy(r, CURVE_Order_BN462);
    BIG_BN462_randtrunc(x1, r, 2 * CURVE_SECURITY_BN462, RNG);
    BIG_BN462_randtrunc(x2, r, 2 * CURVE_SECURITY_BN462, RNG);

    for (i = 0; i < len_sk; ++i){
        BIG_BN462_randtrunc(a[i], r, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(b[i], r, 2 * CURVE_SECURITY_BN462, RNG);
    }

    for (i = 0; i < 2; ++i){
        for (j = 0; j < 2; ++j){
            BIG_BN462_randtrunc(r_u[i][j], r, 2 * CURVE_SECURITY_BN462, RNG);
            BIG_BN462_randtrunc(r_v[i][j], r, 2 * CURVE_SECURITY_BN462, RNG);
            BIG_BN462_randtrunc(r_u_prime[i][j], r, 2 * CURVE_SECURITY_BN462, RNG);
            ECP_BN462_copy(&crs_u[i][j], &gen1);
            ECP2_BN462_copy(&crs_v[i][j], &gen2); 
            ECP_BN462_copy(&crs_u_prime[i][j], &gen1);
        }
    }

    for (i = 0; i < 2; ++i){
        for (j = 0; j < 2; ++j){
            PAIR_BN462_G1mul(&crs_u[i][j], r_u[i][j]);
            PAIR_BN462_G2mul(&crs_v[i][j], r_v[i][j]); 
            PAIR_BN462_G1mul(&crs_u_prime[i][j], r_u[i][j]);
        }
    }

    ECP_BN462_copy(&g, &crs_u[0][0]); 
    ECP_BN462_copy(&h, &crs_u[0][1]);
    ECP_BN462_copy(&g1, &crs_u[1][0]); 
    ECP_BN462_copy(&h1, &crs_u[1][1]);
    ECP2_BN462_copy(&h_hat, &gen2);
    ECP2_BN462_copy(&g_hat, &gen2); 
    PAIR_BN462_G2mul(&h_hat, x1); 
    PAIR_BN462_G2mul(&g_hat, x2);  
   
    ECP_BN462_copy(&f, &g);  
    PAIR_BN462_G1mul(&f, a[0]); 
    ECP_BN462_copy(&temp1, &h); 
    PAIR_BN462_G1mul(&temp1, b[0]); 
    ECP_BN462_add(&f, &temp1); 

    ECP2_BN462_copy(&f1_hat, &g_hat);  
    PAIR_BN462_G2mul(&f1_hat, a[1]); 
    ECP2_BN462_copy(&temp2, &h_hat); 
    PAIR_BN462_G2mul(&temp2, b[1]); 
    ECP2_BN462_add(&f1_hat, &temp2); 
    ECP2_BN462_copy(&f2_hat, &g_hat);  
    PAIR_BN462_G2mul(&f2_hat, a[2]); 
    ECP2_BN462_copy(&temp2, &h_hat); 
    PAIR_BN462_G2mul(&temp2, b[2]); 
    ECP2_BN462_add(&f2_hat, &temp2); 

    for (i = 0; i < len_sig; ++i){
        BIG_BN462_randtrunc(chi_i[i], r, 2 * CURVE_SECURITY_BN462, RNG); 
        BIG_BN462_randtrunc(gamma_i[i], r, 2 * CURVE_SECURITY_BN462, RNG); 
        ECP_BN462_copy(&pk_sxdh[i], &g); 
        PAIR_BN462_G1mul(&pk_sxdh[i], chi_i[i]); 
        ECP_BN462_copy(&temp_i[i], &h); 
        PAIR_BN462_G1mul(&temp_i[i], gamma_i[i]); 
        ECP_BN462_add(&pk_sxdh[i], &temp_i[i]); 
    }

    // Set values to SK
    for (i = 0; i < len_sk; ++i){
        BIG_BN462_rcopy(SK->a[i], a[i]);
        BIG_BN462_rcopy(SK->b[i], b[i]);
    }

    // Set values to PK
    BIG_BN462_rcopy(PK->order, r);
    ECP_BN462_copy(&PK->g, &g);
    ECP_BN462_copy(&PK->h, &h);
    ECP_BN462_copy(&PK->g1, &g1);
    ECP_BN462_copy(&PK->h1, &h1);
    ECP_BN462_copy(&PK->f, &f);

    for (i = 0; i < len_sig; ++i){
        ECP_BN462_copy(&PK->k_i[i], &pk_sxdh[i]);
    }

    for (i = 0; i < 2; ++i){
        for (j = 0; j < 2; ++j){
            ECP_BN462_copy(&PK->crs_u[i][j], &crs_u[i][j]);
            ECP2_BN462_copy(&PK->crs_v[i][j], &crs_v[i][j]);
            ECP_BN462_copy(&PK->crs_u_prime[i][j], &crs_u_prime[i][j]);
        }
    }
    ECP2_BN462_copy(&PK->g_hat, &g_hat);
    ECP2_BN462_copy(&PK->h_hat, &h_hat);
    ECP2_BN462_copy(&PK->f1_hat, &f1_hat);
    ECP2_BN462_copy(&PK->f2_hat, &f2_hat);
}

void encrypt(csprng *RNG, ECP_BN462 M, pk_complex_ballot PK, complex_ballot *ballot){
    int i = 0, j = 0, lenA = 0, lenB = 0, lenC = 0, lenD = 0, lenX = 0, lenY = 0;
    ECP2_BN462 R_hat, S_hat, crhat_1, crhat_2, crhat_0, cshat, temp2, C_R_hat[2], C_S_hat[2];
    ECP_BN462  d1, d2, temp, cm_1, cm_2, cm_0,  C_M[2];
    BIG_BN462 s,r, theta, gamma, R_M[2], R_R_hat[2], R_S_hat[2];

    BIG_BN462_randtrunc(s, PK.order, 2 * CURVE_SECURITY_BN462, RNG);     
    BIG_BN462_randtrunc(r, PK.order, 2 * CURVE_SECURITY_BN462, RNG);
    BIG_BN462_randtrunc(theta, PK.order, 2 * CURVE_SECURITY_BN462, RNG);
    BIG_BN462_randtrunc(gamma, PK.order, 2 * CURVE_SECURITY_BN462, RNG);

    ECP2_BN462_copy(&R_hat, &PK.g_hat); 
    ECP2_BN462_copy(&S_hat, &PK.g_hat); 
    PAIR_BN462_G2mul(&R_hat, r); 
    PAIR_BN462_G2mul(&S_hat, s); 
    ECP_BN462_copy(&d1, &PK.g); 
    PAIR_BN462_G1mul(&d1, r); 
    ECP_BN462_copy(&temp, &PK.h);
    PAIR_BN462_G1mul(&temp, s); 
    ECP_BN462_add(&d1, &temp); 
    ECP_BN462_add(&d1, &M); 
    ECP_BN462_copy(&d2, &PK.g1); 
    PAIR_BN462_G1mul(&d2, r); 
    ECP_BN462_copy(&temp, &PK.h1);
    PAIR_BN462_G1mul(&temp, s); 
    ECP_BN462_add(&d2, &temp); 
    ECP_BN462_copy(&cm_0, &PK.f); 
    PAIR_BN462_G1mul(&cm_0, theta);
    ECP_BN462_add(&cm_0, &M);
    ECP_BN462_copy(&cm_1, &PK.g); 
    PAIR_BN462_G1mul(&cm_1, theta);
    ECP_BN462_copy(&cm_2, &PK.h); 
    PAIR_BN462_G1mul(&cm_2, theta);
    ECP2_BN462_copy(&crhat_0, &PK.f1_hat); 
    PAIR_BN462_G2mul(&crhat_0, gamma);
    ECP2_BN462_add(&crhat_0, &R_hat); 
    ECP2_BN462_copy(&crhat_1, &PK.g_hat); 
    PAIR_BN462_G2mul(&crhat_1, gamma);
    ECP2_BN462_copy(&crhat_2, &PK.h_hat); 
    PAIR_BN462_G2mul(&crhat_2, gamma);
    ECP2_BN462_copy(&cshat, &PK.f2_hat); 
    PAIR_BN462_G2mul(&cshat, gamma);
    ECP2_BN462_add(&cshat, &S_hat);  
  
    for (i = 0; i < 2; ++i){
        BIG_BN462_randtrunc(R_M[i], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(R_R_hat[i], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(R_S_hat[i], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
        
    }
    commit_G1( C_M,  R_M, M,  PK.crs_u);
    commit_G2( C_R_hat,  R_R_hat, R_hat,  PK.crs_v);
    commit_G2( C_S_hat,  R_S_hat, S_hat,  PK.crs_v);

// 1. PROVE THAT CT IS WELL-FORMED: EQUATIONS 10A, 10B
    ECP_BN462 g_inv, h_inv, C_G[2], G;
    ECP2_BN462 C_G_hat[2], G_hat, G_theta, C_theta[2]; 
    BIG_BN462  R_G[2], R_G_hat[2], R_theta[2], S5a[2][2];

    ECP_BN462_copy(&G, &PK.g);
    ECP2_BN462_copy(&G_hat, &PK.g_hat);
    for (i = 0; i < 2; ++i){
        BIG_BN462_randtrunc(R_G[i], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(R_G_hat[i], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(R_theta[i], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
    }
    ECP2_BN462_copy(&G_theta, &G_hat); 
    PAIR_BN462_G2mul(&G_theta, theta); 
    commit_G1( C_G,  R_G, G,  PK.crs_u);
    commit_G2( C_G_hat,  R_G_hat, G_hat,  PK.crs_v);
    commit_G2( C_theta,  R_theta, G_theta,  PK.crs_v);
    ECP_BN462_copy(&g_inv, &PK.g);
    ECP_BN462_neg(&g_inv); 
    ECP_BN462_copy(&h_inv, &PK.h);
    ECP_BN462_neg(&h_inv); 
    ECP_BN462 A5a[2] ={cm_1, g_inv}, pi_eq5a[2], A5b[2] ={cm_2, h_inv}, pi_eq5b[2];
    for (i = 0; i < 2; ++i){
        BIG_464_60_rcopy(S5a[0][i], R_G_hat[i]);
        BIG_464_60_rcopy(S5a[1][i], R_theta[i]);
    }
    lenA = 2;
    prove_pairing_linear_AY(lenA, pi_eq5a, A5a, S5a); 
    prove_pairing_linear_AY(lenA, pi_eq5b, A5b, S5a); 

// 2. PROVE THAT CT IS WELL-FORMED: EQUATION 10E
    ECP_BN462 X[1], qu_temp, f_inv, A[2], theta_eq5e[2][2];
    ECP2_BN462 Y[2], B[1], pi_eq5e[2][2];
    BIG_BN462 Gamma[1][2], R_f[2], R_5[1][2], T[2][2];

    ECP_BN462_copy(&qu_temp, &PK.f); 
    PAIR_BN462_G1mul(&qu_temp, theta);
    ECP_BN462_copy(&X[0], &qu_temp); 
    ECP2_BN462_copy(&Y[0], &G_hat); 
    ECP2_BN462_copy(&Y[1], &G_theta);
    ECP_BN462_copy(&f_inv, &PK.f);
    ECP_BN462_neg(&f_inv);
    iota1(f_inv, A); 
    ECP2_BN462_inf(&B[0]); 
    BIG_464_60_one(Gamma[0][0]); 
    BIG_464_60_zero(Gamma[0][1]);
                                  
    for(i =0; i<2; i++){ 
        BIG_BN462_modneg(R_f[i], R_M[i], PK.order); 
        BIG_464_60_rcopy(R_5[0][i], R_f[i]); 
        for(j =0; j<2; j++){
            BIG_BN462_randtrunc(T[i][j], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
        }
    }

    lenA = 2, lenB = 1, lenX = 1, lenY = 2; 
    prove_pairing_quadratic(lenA, lenB, lenX, lenY, pi_eq5e, theta_eq5e, PK.crs_u, PK.crs_v, X, Y, A, B, Gamma, R_5, S5a, T); 

//3. PROVE THAT CT IS WELL-FORMED: EQUATIONS 10C 10D
    ECP2_BN462 g_hat_inv, pi_11c[2], h_hat_inv, pi_11d[2];;
    ECP_BN462 G_gamma,  C_G_gamma[2];
    BIG_BN462 R_G_gamma[2], R_11c[2][2];

    ECP2_BN462_copy(&g_hat_inv, &PK.g_hat);
    ECP2_BN462_neg(&g_hat_inv);
    ECP2_BN462 B_11c[2] ={crhat_1, g_hat_inv};
    ECP2_BN462_copy(&h_hat_inv, &PK.h_hat);
    ECP2_BN462_neg(&h_hat_inv);
    ECP2_BN462 B_11d[2] ={crhat_2, h_hat_inv};
    ECP_BN462_copy(&G_gamma, &G);
    PAIR_BN462_G1mul(&G_gamma, gamma);

    for (i = 0; i < 2; ++i){
        BIG_BN462_randtrunc(R_G_gamma[i], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
    }

    commit_G1( C_G_gamma,  R_G_gamma, G_gamma,  PK.crs_u);
    for(i =0; i<2; i++){
        BIG_464_60_rcopy(R_11c[0][i], R_G[i]);
        BIG_464_60_rcopy(R_11c[1][i], R_G_gamma[i]);
    }

    lenC = 2;
    prove_pairing_linear_XB(lenC, pi_11c, B_11c, R_11c);
    prove_pairing_linear_XB(lenC, pi_11d, B_11d, R_11c);

//5. PROVE THAT G AND G_hat ARE WELL-FORMED: EQUATION 9
    ECP2_BN462 Y_eq4[1], pi_eq4[2][2], D_4[1][2];
    ECP_BN462 A_eq4[1], theta_eq4[2][2];
    BIG_BN462 Gamma_eq4[1][2], R_4[1][2], S_4[1][2];

    ECP_BN462_copy(&X[0], &G); 
    ECP2_BN462_copy(&Y_eq4[0], &G_hat); 
    ECP_BN462_copy(&A_eq4[0], &PK.g); 
    ECP2_BN462_copy(&B[0], &PK.g_hat); 
    ECP2_BN462_neg(&B[0]);
    
    for(i =0; i<2; i++){ 
        BIG_464_60_zero(Gamma_eq4[0][i]);
        BIG_464_60_rcopy(R_4[0][i], R_G[i]);
         BIG_464_60_rcopy(S_4[0][i], R_G_hat[i]); 
        for(j =0; j<2; j++){
            BIG_BN462_randtrunc(T[i][j], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
        }
    }

    lenA = 1, lenB = 1, lenX = 1, lenY = 1; 
    prove_pairing_quadratic(lenA, lenB, lenX, lenY, pi_eq4, theta_eq4, PK.crs_u, PK.crs_v, X, Y_eq4, A_eq4, B, Gamma_eq4, R_4, S_4, T); 

//6. COMPUTE THE PROOF OF OPENINGS: EQUATION 11A
    ECP2_BN462 Y_12a[3] = {G_hat, S_hat, R_hat}, pi_12a[2][2], B_12a[1];
    ECP_BN462 A_12a[3], theta_12a[2][2], X_12a[1];;
    BIG_BN462 Gamma_12a[1][3], R_12a[1][2], S_12a[3][2];  

    ECP_BN462_copy(&X_12a[0], &M); 
    ECP2_BN462_copy(&B_12a[0], &PK.g_hat); 
    ECP2_BN462_neg(&B_12a[0]);
    ECP_BN462_copy(&A_12a[0], &d1); 
    ECP_BN462_copy(&A_12a[1], &h_inv); 
    ECP_BN462_copy(&A_12a[2], &g_inv); 
    for(i =0; i<3; i++){
        BIG_464_60_zero(Gamma_12a[0][i]);
    }
    
    for(i =0; i<2; i++){ 
        BIG_464_60_rcopy(R_12a[0][i], R_M[i]);
        BIG_464_60_rcopy(S_12a[0][i], R_G_hat[i]); 
        BIG_464_60_rcopy(S_12a[1][i], R_S_hat[i]);
        BIG_464_60_rcopy(S_12a[2][i], R_R_hat[i]); 
        for(j =0; j<2; j++){
            BIG_BN462_randtrunc(T[i][j], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
        }
    }

    lenA = 3, lenB = 1, lenX = 1, lenY = 3;
    prove_pairing_quadratic(lenA, lenB, lenX, lenY, pi_12a, theta_12a, PK.crs_u, PK.crs_v, X_12a, Y_12a, A_12a, B_12a, Gamma_12a, R_12a, S_12a, T); 
    
//7. COMPUTE THE PROOF OF OPENINGS: EQUATION 11B
    ECP_BN462 d2_inv;
    BIG_BN462 S12b[3][2];
    ECP_BN462_copy(&d2_inv, &d2);
    ECP_BN462_neg(&d2_inv);
    ECP_BN462 A12b[3] ={PK.g1, PK.h1, d2_inv}, pi_12b[2];

    for (i = 0; i < 2; ++i){
        BIG_464_60_rcopy(S12b[0][i], R_R_hat[i]);
        BIG_464_60_rcopy(S12b[1][i], R_S_hat[i]);
        BIG_464_60_rcopy(S12b[2][i], R_G_hat[i]);
    }
    
    lenA = 3;
    prove_pairing_linear_AY(lenA, pi_12b, A12b, S12b); 
 
//8. COMPUTE THE PROOF OF TRACEABILITY: EQUATION 13
    int slen=3;
    BIG_BN462 chi_trace[slen], gamma_trace[slen], R_Z1[2], R_R1[2], R_14[2][2];
    ECP_BN462  Z_sig1, R_sig1, Z1[1], R1[1], Z2[1], R2[1], Z3[1], R3[1], msg1[slen], msg2[slen], msg3[slen], C_Z1[2], C_R1[2];
    ECP2_BN462 pk_sig[slen], temp2_i[slen], pi_14[2], B_14[2] ={PK.g_hat, PK.h_hat};
    
    for (i = 0; i < slen; ++i){
        BIG_BN462_randtrunc(chi_trace[i], PK.order, 2 * CURVE_SECURITY_BN462, RNG); 
        BIG_BN462_randtrunc(gamma_trace[i], PK.order, 2 * CURVE_SECURITY_BN462, RNG); 
        ECP2_BN462_copy(&pk_sig[i], &PK.g_hat); 
        ECP2_BN462_copy(&temp2_i[i], &PK.h_hat); 
        PAIR_BN462_G2mul(&temp2_i[i], gamma_trace[i]); 
        PAIR_BN462_G2mul(&pk_sig[i], chi_trace[i]); 
        ECP2_BN462_add(&pk_sig[i], &temp2_i[i]); 
    }
    // sign row 1
    ECP_BN462_copy(&msg1[0], &PK.g);
    ECP_BN462_copy(&msg1[1], &d1);
    ECP_BN462_copy(&msg1[2], &d2);
    ECP_BN462_inf(&Z1[0]);
    ECP_BN462_inf(&R1[0]);
    sign_lhspG1(slen, Z1,  R1,  gamma_trace,  chi_trace, msg1);   

    // sign row 2 
    ECP_BN462_inf(&msg2[0]);
    ECP_BN462_copy(&msg2[1], &PK.g);
    ECP_BN462_copy(&msg2[2], &PK.g1);
    ECP_BN462_inf(&Z2[0]);
    ECP_BN462_inf(&R2[0]);
    sign_lhspG1(slen, Z2,  R2,  gamma_trace,  chi_trace, msg2);

    // sign row 3 
    ECP_BN462_inf(&msg3[0]);
    ECP_BN462_copy(&msg3[1], &PK.h);
    ECP_BN462_copy(&msg3[2], &PK.h1);
    ECP_BN462_inf(&Z3[0]);
    ECP_BN462_inf(&R3[0]);
    sign_lhspG1(slen, Z3,  R3,  gamma_trace,  chi_trace, msg3);

    for (i = 0; i < 2; ++i){
        BIG_BN462_randtrunc(R_Z1[i], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(R_R1[i], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
    }

    ECP_BN462_copy(&Z_sig1, &Z1[0]);
    ECP_BN462_copy(&R_sig1, &R1[0]);
    commit_G1( C_Z1,  R_Z1, Z_sig1,  PK.crs_u_prime); 
    commit_G1( C_R1,  R_R1, R_sig1,  PK.crs_u_prime);
    
    for (i = 0; i < 2; ++i){
        BIG_464_60_rcopy(R_14[0][i], R_Z1[i]);
        BIG_464_60_rcopy(R_14[1][i], R_R1[i]);
    }
    
    lenC = 2;
    prove_pairing_linear_XB(lenC, pi_14, B_14, R_14);

//9.  COMPUTE THE TCCA PROOF: EQUATION 14
    ECP2_BN462 tau, C_X[2], A_inf, B_inf, C_A[2], C_B[2];
    ECP_BN462 pk_sxdh, pi_ss[2], A_ss[3] ={PK.g, PK.h, pk_sxdh};; 
    BIG_BN462 R_A[2], R_B[2], R_X[2], tau_big,  R_ss[3][2];
    char oct[128];
    octet  tau_oct = {1, sizeof(PK.order), oct}; 

    ECP2_BN462_inf(&tau);
    ECP2_BN462_inf(&A_inf);
    ECP2_BN462_inf(&B_inf);
    derive_commitment2(C_X,  C_G_hat,  PK.g_hat); 
   
    for (i = 0; i < 2; ++i){
        BIG_BN462_modneg(R_X[i], R_G_hat[i], PK.order);
        BIG_BN462_randtrunc(R_A[i], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(R_B[i], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
    }
    commit_G2( C_A,  R_A, A_inf,  PK.crs_v);
    commit_G2( C_B,  R_B, B_inf,  PK.crs_v);
    
    for (i = 0; i < 3; ++i){
        ECP2_BN462_add(&tau, &pk_sig[i]);
    }
    ECP2_BN462_toOctet(&tau_oct, &tau, true);
    BIG_464_60_fromBytes(tau_big, tau_oct.val);
    ECP_BN462_copy(&pk_sxdh, &PK.k_i[1]);
    PAIR_BN462_G1mul(&pk_sxdh, tau_big);
    ECP_BN462_add(&pk_sxdh, &PK.k_i[0]);
    ECP_BN462_neg(&pk_sxdh);
    ECP_BN462_copy(&A_ss[2], &pk_sxdh);
    
    for (i = 0; i < 2; ++i){
        BIG_464_60_rcopy(R_ss[0][i], R_A[i]);
        BIG_464_60_rcopy(R_ss[1][i], R_B[i]);
        BIG_464_60_rcopy(R_ss[2][i], R_X[i]);
    }
    
    lenA = 3;
    prove_pairing_linear_AY(lenA, pi_ss, A_ss, R_ss); 
    
//10. PROVE THAT CT IS WELL-FORMED: EQUATION 10F
    ECP_BN462 X_11f[2], A_11f[1], theta_11f[2][2];
    ECP2_BN462 Y_11f[1], B_11f[2], pi_11f[2][2], qu_temp2, f_hat_inv;
    BIG_BN462 Gamma_11f[2][1], R_f1_hat[2], R_11f[2][2], S_11f[1][2]; 

    ECP_BN462_copy(&X_11f[0], &G); 
    ECP_BN462_copy(&X_11f[1], &G_gamma); 
    ECP2_BN462_copy(&qu_temp2, &PK.f1_hat); 
    PAIR_BN462_G2mul(&qu_temp2, gamma);
    ECP2_BN462_copy(&Y_11f[0], &qu_temp2);
    ECP2_BN462_copy(&f_hat_inv, &PK.f1_hat);
    ECP2_BN462_neg(&f_hat_inv);
    iota2(f_hat_inv, B_11f); 
    ECP_BN462_inf(&A_11f[0]); 
    BIG_464_60_one(Gamma_11f[0][0]); 
    BIG_464_60_zero(Gamma_11f[1][0]);
                            
    for(i =0; i<2; i++){ 
        BIG_BN462_modneg(R_f1_hat[i], R_R_hat[i], PK.order); 
        BIG_464_60_rcopy(R_11f[0][i], R_G[i]); 
        BIG_464_60_rcopy(R_11f[1][i], R_G_gamma[i]); 
        BIG_464_60_rcopy(S_11f[0][i], R_f1_hat[i]);
        for(j =0; j<2; j++){
            BIG_BN462_randtrunc(T[i][j], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
        }
    }

    lenA = 1, lenB = 2, lenX = 2, lenY = 1; 
    prove_pairing_quadratic2(lenA, lenB, lenX, lenY, pi_11f, theta_11f, PK.crs_u, PK.crs_v, X_11f, Y_11f, A_11f, B_11f, Gamma_11f, R_11f, S_11f, T); 

//10. PROVE THAT CT IS WELL-FORMED: EQUATION 10G
    ECP_BN462 X_11g[2], A_11g[1], theta_11g[2][2];
    ECP2_BN462 Y_11g[1], B_11g[2], pi_11g[2][2];
    BIG_BN462 Gamma_11g[2][1], R_f2_hat[2], R_11g[2][2], S_11g[1][2]; 

    ECP_BN462_copy(&X_11g[0], &G); 
    ECP_BN462_copy(&X_11g[1], &G_gamma); 
    ECP2_BN462_copy(&qu_temp2, &PK.f2_hat); 
    PAIR_BN462_G2mul(&qu_temp2, gamma);
    ECP2_BN462_copy(&Y_11g[0], &qu_temp2);
    ECP2_BN462_copy(&f_hat_inv, &PK.f2_hat);
    ECP2_BN462_neg(&f_hat_inv);
    iota2(f_hat_inv, B_11g); 
    ECP_BN462_inf(&A_11g[0]); 
    BIG_464_60_one(Gamma_11g[0][0]); 
    BIG_464_60_zero(Gamma_11g[1][0]);
                             
    for(i =0; i<2; i++){ 
        BIG_BN462_modneg(R_f2_hat[i], R_S_hat[i], PK.order);  
        BIG_464_60_rcopy(R_11g[0][i], R_G[i]); 
        BIG_464_60_rcopy(R_11g[1][i], R_G_gamma[i]);
        BIG_464_60_rcopy(S_11g[0][i], R_f2_hat[i]);  
        for(j =0; j<2; j++){
            BIG_BN462_randtrunc(T[i][j], PK.order, 2 * CURVE_SECURITY_BN462, RNG);
        }
    }

    lenA = 1, lenB = 2, lenX = 2, lenY = 1;
    prove_pairing_quadratic2(lenA, lenB, lenX, lenY, pi_11g, theta_11g, PK.crs_u, PK.crs_v, X_11g, Y_11g, A_11g, B_11g, Gamma_11g, R_11g, S_11g, T); 

    // Set values to ballot
    ECP_BN462_copy(&ballot->cm_0, &cm_0);
    ECP_BN462_copy(&ballot->cm_1, &cm_1);
    ECP_BN462_copy(&ballot->cm_2, &cm_2);
    ECP_BN462_copy(&ballot->d1, &d1);
    ECP_BN462_copy(&ballot->d2, &d2);

    ECP2_BN462_copy(&ballot->crhat_0, &crhat_0);
    ECP2_BN462_copy(&ballot->crhat_1, &crhat_1);
    ECP2_BN462_copy(&ballot->crhat_2, &crhat_2);
    ECP2_BN462_copy(&ballot->cshat, &cshat);
    for(i =0; i<slen; i++){
        ECP2_BN462_copy(&ballot->pk_sig[i], &pk_sig[i]);
    }
    for(i =0; i<2; i++){
        ECP_BN462_copy(&ballot->C_M[i], &C_M[i]);
        ECP_BN462_copy(&ballot->C_G[i], &C_G[i]);
        ECP_BN462_copy(&ballot->C_G_gamma[i], &C_G_gamma[i]);
        ECP_BN462_copy(&ballot->C_Z1[i], &C_Z1[i]);
        ECP_BN462_copy(&ballot->C_R1[i], &C_R1[i]);
        ECP2_BN462_copy(&ballot->C_R_hat[i], &C_R_hat[i]);
        ECP2_BN462_copy(&ballot->C_S_hat[i], &C_S_hat[i]);
        ECP2_BN462_copy(&ballot->C_G_hat[i], &C_G_hat[i]);
        ECP2_BN462_copy(&ballot->C_theta[i], &C_theta[i]);
        ECP2_BN462_copy(&ballot->C_A[i], &C_A[i]);
        ECP2_BN462_copy(&ballot->C_B[i], &C_B[i]);
        ECP2_BN462_copy(&ballot->pi_11c[i], &pi_11c[i]);
        ECP2_BN462_copy(&ballot->pi_11d[i], &pi_11d[i]);
        ECP2_BN462_copy(&ballot->pi_14[i], &pi_14[i]);
        ECP_BN462_copy(&ballot->pi_eq5a[i], &pi_eq5a[i]);
        ECP_BN462_copy(&ballot->pi_eq5b[i], &pi_eq5b[i]); 
        ECP_BN462_copy(&ballot->pi_12b[i], &pi_12b[i]); 
        ECP_BN462_copy(&ballot->pi_ss[i], &pi_ss[i]); 
        for(j =0; j<2; j++){
            ECP2_BN462_copy(&ballot->pi_eq5e[i][j], &pi_eq5e[i][j]);
            ECP_BN462_copy(&ballot->theta_eq5e[i][j], &theta_eq5e[i][j]);
            ECP2_BN462_copy(&ballot->pi_eq4[i][j], &pi_eq4[i][j]);
            ECP_BN462_copy(&ballot->theta_eq4[i][j], &theta_eq4[i][j]);
            ECP2_BN462_copy(&ballot->pi_12a[i][j], &pi_12a[i][j]);
            ECP_BN462_copy(&ballot->theta_12a[i][j], &theta_12a[i][j]); 
            ECP2_BN462_copy(&ballot->pi_11f[i][j], &pi_11f[i][j]);
            ECP_BN462_copy(&ballot->theta_11f[i][j], &theta_11f[i][j]);
            ECP2_BN462_copy(&ballot->pi_11g[i][j], &pi_11g[i][j]);
            ECP_BN462_copy(&ballot->theta_11g[i][j], &theta_11g[i][j]);
        } 
    }
}

bool verify(pk_complex_ballot PK, complex_ballot ballot){
    int i = 0, lenA = 0, lenB = 0, lenC = 0, lenD = 0, lenX = 0, lenY = 0;

//1. 2. VERIFY EQUATIONS 10A 10B
    ECP2_BN462 D5a[2][2];
    ECP_BN462 g_inv, h_inv;
    FP12_BN462 t_eq5a, t_eq5b;

    ECP_BN462_copy(&g_inv, &PK.g);
    ECP_BN462_neg(&g_inv); 
    ECP_BN462_copy(&h_inv, &PK.h);
    ECP_BN462_neg(&h_inv); 
    ECP_BN462 A5a[2] ={ballot.cm_1, g_inv}, A5b[2] ={ballot.cm_2, h_inv};
    lenA = 2;

    for(i =0; i<2; i++){
        ECP2_BN462_copy(&D5a[0][i], &ballot.C_G_hat[i]);
        ECP2_BN462_copy(&D5a[1][i], &ballot.C_theta[i]);
    }

    FP12_BN462_one(&t_eq5a); 
    if(!verify_pairing_linear_AY(lenA, PK.crs_v, A5a, D5a, ballot.pi_eq5a, t_eq5a))
        return false;

    FP12_BN462_one(&t_eq5b);
    if(!verify_pairing_linear_AY(lenA, PK.crs_v, A5b, D5a, ballot.pi_eq5b, t_eq5b))
        return false;

//3. VERIFY EQUATION 10E
    FP12_BN462 t_eq5e;
    ECP_BN462  C_f[2], C[1][2], A[2], f_inv;
    ECP2_BN462 D[2][2], B[1];
    BIG_BN462 Gamma[1][2];

    ECP_BN462_copy(&f_inv, &PK.f);
    ECP_BN462_neg(&f_inv);
    iota1(f_inv, A); 
    ECP2_BN462_inf(&B[0]); 
    BIG_464_60_one(Gamma[0][0]); 
    BIG_464_60_zero(Gamma[0][1]);
    FP12_BN462_one(&t_eq5e);
    derive_commitment1(C_f,  ballot.C_M,  ballot.cm_0);

    for(i =0; i<2; i++){
        ECP_BN462_copy(&C[0][i], &C_f[i]);
        ECP2_BN462_copy(&D[0][i], &ballot.C_G_hat[i]);
        ECP2_BN462_copy(&D[1][i], &ballot.C_theta[i]);
    }
    
    lenA =2, lenC = 1, lenD = 2;
    if(!verify_pairing_quadratic1(lenA, lenC, lenD, PK.crs_u, PK.crs_v, C, D, A, B, Gamma, ballot.pi_eq5e, ballot.theta_eq5e, t_eq5e))
        return false; 

//4.5. VERIFY EQUATIONS 10C 10D
    FP12_BN462 t_11c, t_11d;
    ECP_BN462 C_11c[2][2];
    ECP2_BN462 g_hat_inv, h_hat_inv;

    ECP2_BN462_copy(&g_hat_inv, &PK.g_hat);
    ECP2_BN462_neg(&g_hat_inv);
    ECP2_BN462 B_11c[2] ={ballot.crhat_1, g_hat_inv};
    ECP2_BN462_copy(&h_hat_inv, &PK.h_hat);
    ECP2_BN462_neg(&h_hat_inv);
    ECP2_BN462 B_11d[2] ={ballot.crhat_2, h_hat_inv};

    for(i =0; i<2; i++){
        ECP_BN462_copy(&C_11c[0][i], &ballot.C_G[i]);
        ECP_BN462_copy(&C_11c[1][i], &ballot.C_G_gamma[i]);
    }

    lenC = 2;
    FP12_BN462_one(&t_11c);
    if(! verify_pairing_linear_XB(lenC, PK.crs_u, B_11c, C_11c, ballot.pi_11c, t_11c))
        return false;

    FP12_BN462_one(&t_11d);
    if(! verify_pairing_linear_XB(lenC, PK.crs_u, B_11d, C_11c, ballot.pi_11d, t_11d))
        return false;

//6. VERIFY EQUATION 9
    FP12_BN462 t_eq4;
    ECP2_BN462 D_4[1][2];
    ECP_BN462 A_eq4[1];
    BIG_BN462 Gamma_eq4[1][2];

    FP12_BN462_one(&t_eq4);
    ECP_BN462_copy(&A_eq4[0], &PK.g); 
    ECP2_BN462_copy(&B[0], &PK.g_hat); 
    ECP2_BN462_neg(&B[0]);
    for(i =0; i<2; i++){
        BIG_464_60_zero(Gamma_eq4[0][i]);
        ECP_BN462_copy(&C[0][i], &ballot.C_G[i]);
        ECP2_BN462_copy(&D_4[0][i], &ballot.C_G_hat[i]);
    }

    lenA =1, lenC = 1, lenD = 1;
    if(!verify_pairing_quadratic(lenA, lenC, lenD, PK.crs_u, PK.crs_v, C, D_4, A_eq4, B, Gamma_eq4, ballot.pi_eq4, ballot.theta_eq4, t_eq4))
        return false;

//7. VERIFY EQUATION 11A
    FP12_BN462 t_12a;
    ECP2_BN462 D_12a[3][2], B_12a[1];
    ECP_BN462 A_12a[3];
    BIG_BN462 Gamma_12a[1][3];

    FP12_BN462_one(&t_12a);
    ECP_BN462_copy(&A_12a[0], &ballot.d1); 
    ECP_BN462_copy(&A_12a[1], &h_inv); 
    ECP_BN462_copy(&A_12a[2], &g_inv);
    ECP2_BN462_copy(&B_12a[0], &PK.g_hat); 
    ECP2_BN462_neg(&B_12a[0]);

    for (i = 0; i < 3; ++i){
        BIG_464_60_zero(Gamma_12a[0][i]);
    }
    for (i = 0; i < 2; ++i){
        ECP_BN462_copy(&C[0][i], &ballot.C_M[i]);
        ECP2_BN462_copy(&D_12a[0][i], &ballot.C_G_hat[i]);
        ECP2_BN462_copy(&D_12a[1][i], &ballot.C_S_hat[i]);
        ECP2_BN462_copy(&D_12a[2][i], &ballot.C_R_hat[i]);
    }
    
    lenA = 3, lenC = 1, lenD = 3;
    if(!verify_pairing_quadratic(lenA, lenC, lenD, PK.crs_u, PK.crs_v, C, D_12a, A_12a, B_12a, Gamma_12a, ballot.pi_12a, ballot.theta_12a, t_12a))
        return false;

//8. VERIFY EQUATION 11B
    ECP2_BN462 D12b[3][2];
    FP12_BN462 t_12b;
    ECP_BN462 d2_inv;

    ECP_BN462_copy(&d2_inv, &ballot.d2);
    ECP_BN462_neg(&d2_inv);
    ECP_BN462 A12b[3] ={PK.g1, PK.h1, d2_inv};
    for (i = 0; i < 2; ++i){
        ECP2_BN462_copy(&D12b[0][i], &ballot.C_R_hat[i]);
        ECP2_BN462_copy(&D12b[1][i], &ballot.C_S_hat[i]);
        ECP2_BN462_copy(&D12b[2][i], &ballot.C_G_hat[i]);
    }
    
    FP12_BN462_one(&t_12b); 
    lenA =3;
    if(! verify_pairing_linear_AY(lenA, PK.crs_v, A12b, D12b, ballot.pi_12b, t_12b))
        return false;

//9. VERIFY EQUATION 13
    ECP2_BN462 B_14[2] ={PK.g_hat, PK.h_hat};
    FP12_BN462 t_14, temp_sig;
    ECP_BN462 C_14[2][2], msg1[3];

    for (i = 0; i < 2; ++i){
        ECP_BN462_copy(&C_14[0][i], &ballot.C_Z1[i]);
        ECP_BN462_copy(&C_14[1][i], &ballot.C_R1[i]);
    }
    FP12_BN462_one(&t_14); 
    ECP_BN462_copy(&msg1[0], &PK.g);
    ECP_BN462_copy(&msg1[1], &ballot.d1);
    ECP_BN462_copy(&msg1[2], &ballot.d2);
    for (i = 0; i < 3; ++i){
        temp_sig= PAIR_BN462_e(ballot.pk_sig[i], msg1[i]);
        FP12_BN462_ssmul(&t_14,&temp_sig);
    }

    lenC = 2;
    if(! verify_pairing_linear_XB(lenC, PK.crs_u_prime, B_14, C_14, ballot.pi_14, t_14))
        return false;

//10. VERIFY EQUATION 14
    ECP2_BN462 D_ss[3][2], tau, C_X[2];
    FP12_BN462 t_ss;
    ECP_BN462 pk_sxdh, A_ss[3] ={PK.g, PK.h, pk_sxdh};;
    BIG_BN462 tau_big;
    char oct[128];
    octet  tau_oct = {1, sizeof(PK.order), oct}; 

    ECP2_BN462_inf(&tau);
    for (i = 0; i < 3; ++i){
        ECP2_BN462_add(&tau, &ballot.pk_sig[i]);
    }
    ECP2_BN462_toOctet(&tau_oct, &tau, true);
    BIG_464_60_fromBytes(tau_big, tau_oct.val);
    ECP_BN462_copy(&pk_sxdh, &PK.k_i[1]);
    PAIR_BN462_G1mul(&pk_sxdh, tau_big);
    ECP_BN462_add(&pk_sxdh, &PK.k_i[0]);
    ECP_BN462_neg(&pk_sxdh);
    ECP_BN462_copy(&A_ss[2], &pk_sxdh);
    derive_commitment2(C_X,  ballot.C_G_hat,  PK.g_hat); 

    for (i = 0; i < 2; ++i){
        ECP2_BN462_copy(&D_ss[0][i], &ballot.C_A[i]);
        ECP2_BN462_copy(&D_ss[1][i], &ballot.C_B[i]);
        ECP2_BN462_copy(&D_ss[2][i], &C_X[i]);
    }
   
    FP12_BN462_one(&t_ss);
    lenA = 3; 
    if(! verify_pairing_linear_AY(lenA, PK.crs_v, A_ss, D_ss, ballot.pi_ss, t_ss))
        return false;

//11. VERIFY EQUATION 10F
    FP12_BN462 t_11f;
    ECP_BN462  C_11f[2][2], A_11f[1];
    ECP2_BN462 D_11f[1][2], D_f1_hat[2], B_11f[2], f_hat_inv;
    BIG_BN462 Gamma_11f[2][1];  

    ECP_BN462_inf(&A_11f[0]); 
    ECP2_BN462_copy(&f_hat_inv, &PK.f1_hat);
    ECP2_BN462_neg(&f_hat_inv);
    iota2(f_hat_inv, B_11f);
    BIG_464_60_one(Gamma_11f[0][0]); 
    BIG_464_60_zero(Gamma_11f[1][0]); 
    FP12_BN462_one(&t_11f);
    derive_commitment2(D_f1_hat,  ballot.C_R_hat,  ballot.crhat_0); 

    for (i = 0; i < 2; ++i){
        ECP_BN462_copy(&C_11f[0][i], &ballot.C_G[i]);
        ECP_BN462_copy(&C_11f[1][i], &ballot.C_G_gamma[i]);
        ECP2_BN462_copy(&D_11f[0][i], &D_f1_hat[i]);
    }
    lenA = 1, lenC = 2, lenD = 1;
    if(! verify_pairing_quadratic2(lenA, lenC, lenD, PK.crs_u, PK.crs_v, C_11f, D_11f, A_11f, B_11f, Gamma_11f, ballot.pi_11f, ballot.theta_11f, t_11f))
        return false;

//12. VERIFY EQUATION 10G
    FP12_BN462 t_11g;
    ECP_BN462  C_11g[2][2], A_11g[1];
    ECP2_BN462 D_11g[1][2], D_f2_hat[2], B_11g[2];
    BIG_BN462 Gamma_11g[2][1];

    derive_commitment2(D_f2_hat,  ballot.C_S_hat,  ballot.cshat); 
    ECP_BN462_inf(&A_11g[0]); 
    ECP2_BN462_copy(&f_hat_inv, &PK.f2_hat);
    ECP2_BN462_neg(&f_hat_inv);
    iota2(f_hat_inv, B_11g);
    BIG_464_60_one(Gamma_11g[0][0]); 
    BIG_464_60_zero(Gamma_11g[1][0]);
    FP12_BN462_one(&t_11g);
    derive_commitment2(D_f2_hat,  ballot.C_S_hat,  ballot.cshat); 

    for (i = 0; i < 2; ++i){
        ECP_BN462_copy(&C_11g[0][i], &ballot.C_G[i]);
        ECP_BN462_copy(&C_11g[1][i], &ballot.C_G_gamma[i]);
        ECP2_BN462_copy(&D_11g[0][i], &D_f2_hat[i]);
    }
    lenA = 1, lenC = 2, lenD = 1;
    if(! verify_pairing_quadratic2(lenA, lenC, lenD, PK.crs_u, PK.crs_v, C_11g, D_11g, A_11g, B_11g, Gamma_11g, ballot.pi_11g, ballot.theta_11g, t_11g))
        return false;
    
    return true;
}

void decrypt(ECP_BN462 *decrypted_M, sk_complex_ballot SK, pk_complex_ballot PK, complex_ballot ballot){
    if(! verify(PK, ballot)){
        printf("Verification returns 0!\n");
    } 
    else{
        ECP_BN462 M, temp;
        ECP_BN462_copy(&M,&ballot.cm_0);
        ECP_BN462_copy(&temp,&ballot.cm_1);
        PAIR_BN462_G1mul(&temp,SK.a[0]);
        ECP_BN462_neg(&temp);
        ECP_BN462_add(&M,&temp);
        ECP_BN462_copy(&temp,&ballot.cm_2);
        PAIR_BN462_G1mul(&temp,SK.b[0]);
        ECP_BN462_neg(&temp);
        ECP_BN462_add(&M,&temp);
        ECP_BN462_copy(decrypted_M, &M);
    }
}



