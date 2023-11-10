#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "tools.h"

void keygen(csprng *RNG, sk_simple_ballot *SK, pk_simple_ballot *PK){
    int i = 0, j = 0, len_sk = 3, len_sig = 2;
    ECP_BN462 gen1, g, g1, h1, f1, f2, f3, temp, crs_u[2][2];;
    ECP2_BN462  gen2, g_hat, h_hat, g1_hat, g2_hat, h1_hat, h2_hat, crs_v[2][2], crs_v_prime[2][2], g_hat_i[len_sig], temp_i[len_sig], f_hat[len_sig];
    BIG_BN462 r, a[len_sk], b[len_sk], x[5], r_u[2][2], r_v[2][2], r_v_prime[2][2], chi_i[len_sig], gamma_i[len_sig];

    ECP_BN462_generator(&gen1);
    ECP2_BN462_generator(&gen2);
    BIG_BN462_rcopy(r, CURVE_Order_BN462);
    ECP_BN462_copy(&g1, &gen1); 
    ECP_BN462_copy(&h1, &gen1); 
    ECP_BN462_copy(&g, &gen1); 
    ECP2_BN462_copy(&h_hat, &gen2);
    ECP2_BN462_copy(&g_hat, &gen2);

    for (i = 0; i < len_sk; ++i){
        BIG_BN462_randtrunc(a[i], r, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(b[i], r, 2 * CURVE_SECURITY_BN462, RNG);
    }
    for (i = 0; i < 5; ++i){
        BIG_BN462_randtrunc(x[i], r, 2 * CURVE_SECURITY_BN462, RNG);
    }
   
    for (i = 0; i < 2; ++i){
        for (j = 0; j < 2; ++j){
            BIG_BN462_randtrunc(r_u[i][j], r, 2 * CURVE_SECURITY_BN462, RNG);
            BIG_BN462_randtrunc(r_v[i][j], r, 2 * CURVE_SECURITY_BN462, RNG);
            BIG_BN462_randtrunc(r_v_prime[i][j], r, 2 * CURVE_SECURITY_BN462, RNG);
        }
    }

    for (i = 0; i < 2; ++i){
        for (j = 0; j < 2; ++j){
            ECP_BN462_copy(&crs_u[i][j], &gen1);
            ECP2_BN462_copy(&crs_v[i][j], &gen2);
            ECP2_BN462_copy(&crs_v_prime[i][j], &gen2);
            PAIR_BN462_G1mul(&crs_u[i][j], r_u[i][j]);
            PAIR_BN462_G2mul(&crs_v[i][j], r_v[i][j]); 
            PAIR_BN462_G2mul(&crs_v_prime[i][j], r_v_prime[i][j]); 
        }
    }
    
    PAIR_BN462_G1mul(&g, x[0]); 
    PAIR_BN462_G1mul(&g1, x[1]); 
    PAIR_BN462_G1mul(&h1, x[2]); 
    PAIR_BN462_G2mul(&h_hat, x[3]); 
    PAIR_BN462_G2mul(&g_hat, x[4]); 
    PAIR_BN462_G1mul(&h1, x[2]); 
    ECP_BN462_copy(&f1, &g1);  
    PAIR_BN462_G1mul(&f1, a[0]); 
    ECP_BN462_copy(&temp, &h1); 
    PAIR_BN462_G1mul(&temp, b[0]); 
    ECP_BN462_add(&f1, &temp); 
    ECP_BN462_copy(&f2, &g1);  
    PAIR_BN462_G1mul(&f2, a[1]); 
    ECP_BN462_copy(&temp, &h1); 
    PAIR_BN462_G1mul(&temp, b[1]); 
    ECP_BN462_add(&f2, &temp); 
    ECP_BN462_copy(&f3, &g1);  
    PAIR_BN462_G1mul(&f3, a[2]); 
    ECP_BN462_copy(&temp, &h1); 
    PAIR_BN462_G1mul(&temp, b[2]); 
    ECP_BN462_add(&f3, &temp); 
    
    ECP2_BN462_copy(&g1_hat, &crs_v[0][0]); 
    ECP2_BN462_copy(&h1_hat, &crs_v[0][1]);
    ECP2_BN462_copy(&g2_hat, &crs_v[1][0]); 
    ECP2_BN462_copy(&h2_hat, &crs_v[1][1]);

    for (int i = 0; i < len_sig; ++i){
        BIG_BN462_randtrunc(chi_i[i], r, 2 * CURVE_SECURITY_BN462, RNG); 
        BIG_BN462_randtrunc(gamma_i[i], r, 2 * CURVE_SECURITY_BN462, RNG); 
        ECP2_BN462_copy(&g_hat_i[i], &g_hat); 
        ECP2_BN462_copy(&temp_i[i], &h_hat); 
        PAIR_BN462_G2mul(&temp_i[i], gamma_i[i]); 
        PAIR_BN462_G2mul(&g_hat_i[i], chi_i[i]); 
        ECP2_BN462_add(&g_hat_i[i], &temp_i[i]); 
    }
    for (int i = 0; i < len_sig; ++i){
        ECP2_BN462_copy(&f_hat[i], &g_hat_i[i]); 
    }

    // Set values to SK
    for (int i = 0; i < len_sk; ++i){
        BIG_BN462_rcopy(SK->alpha[i], a[i]);
        BIG_BN462_rcopy(SK->beta[i], b[i]);
    }

    // Set values to PK
    BIG_BN462_rcopy(PK->r, r);
    ECP_BN462_copy(&PK->g, &g);
    ECP_BN462_copy(&PK->g1, &g1);
    ECP_BN462_copy(&PK->h1, &h1);
    ECP_BN462_copy(&PK->f1, &f1);
    ECP_BN462_copy(&PK->f2, &f2);
    ECP_BN462_copy(&PK->f3, &f3);
    for (i = 0; i < 2; ++i){
        ECP2_BN462_copy(&PK->f_hat[i], &f_hat[i]);
        for (j = 0; j < 2; ++j){
            ECP_BN462_copy(&PK->crs_u[i][j], &crs_u[i][j]);
            ECP2_BN462_copy(&PK->crs_v[i][j], &crs_v[i][j]);
            ECP2_BN462_copy(&PK->crs_v_prime[i][j], &crs_v_prime[i][j]);
        }
    }
    ECP2_BN462_copy(&PK->g_hat, &g_hat);
    ECP2_BN462_copy(&PK->h_hat, &h_hat);
    ECP2_BN462_copy(&PK->g1_hat, &g1_hat);
    ECP2_BN462_copy(&PK->g2_hat, &g2_hat);
    ECP2_BN462_copy(&PK->h1_hat, &h1_hat);
    ECP2_BN462_copy(&PK->h2_hat, &h2_hat);
}

void encrypt(csprng *RNG, BIG_BN462 m, pk_simple_ballot PK, simple_ballot *ballot){
    int i = 0, j = 0, lenA = 0, lenB = 0, lenC = 0, lenD = 0, lenX = 0, lenY = 0;
    ECP_BN462 M, R, Q, c_1, c_2, c_0, c_r, c_q, temp1,  C_M[2], C_R[2], C_Q[2], G, C_G[2], g1_inv, h1_inv;
    ECP2_BN462  d1_hat, d2_hat, temp2, H_hat,H_theta, C_H[2], C_theta[2];;
    BIG_BN462 s,q, theta, R_M[2], R_R[2], R_Q[2], R_G[2], R_H[2], R_theta[2];

    BIG_BN462_randtrunc(s, PK.r, 2 * CURVE_SECURITY_BN462, RNG);     
    BIG_BN462_randtrunc(q, PK.r, 2 * CURVE_SECURITY_BN462, RNG);
    BIG_BN462_randtrunc(theta, PK.r, 2 * CURVE_SECURITY_BN462, RNG);

    ECP_BN462_copy(&M, &PK.g); 
    ECP_BN462_copy(&R, &PK.g); 
    ECP_BN462_copy(&Q, &PK.g); 
    PAIR_BN462_G1mul(&M, m); 
    PAIR_BN462_G1mul(&R, s); 
    PAIR_BN462_G1mul(&Q, q); 
    ECP2_BN462_copy(&d1_hat, &PK.g_hat); 
    PAIR_BN462_G2mul(&d1_hat, m); 
    ECP2_BN462_copy(&temp2, &PK.g1_hat);
    PAIR_BN462_G2mul(&temp2, s); 
    ECP2_BN462_add(&d1_hat, &temp2); 
    ECP2_BN462_copy(&temp2, &PK.h1_hat);
    PAIR_BN462_G2mul(&temp2, q); 
    ECP2_BN462_add(&d1_hat, &temp2); 
    ECP2_BN462_copy(&d2_hat, &PK.g2_hat); 
    PAIR_BN462_G2mul(&d2_hat, s); 
    ECP2_BN462_copy(&temp2, &PK.h2_hat);
    PAIR_BN462_G2mul(&temp2, q); 
    ECP2_BN462_add(&d2_hat, &temp2);
    ECP_BN462_copy(&c_1, &PK.g1); 
    PAIR_BN462_G1mul(&c_1, theta);
    ECP_BN462_copy(&c_2, &PK.h1); 
    PAIR_BN462_G1mul(&c_2, theta);
    ECP_BN462_copy(&temp1, &PK.f1);
    PAIR_BN462_G1mul(&temp1, theta); 
    ECP_BN462_copy(&c_0, &M); 
    ECP_BN462_add(&c_0, &temp1); 
    ECP_BN462_copy(&temp1, &PK.f2);
    PAIR_BN462_G1mul(&temp1, theta); 
    ECP_BN462_copy(&c_r, &R); 
    ECP_BN462_add(&c_r, &temp1);
    ECP_BN462_copy(&temp1, &PK.f3);
    PAIR_BN462_G1mul(&temp1, theta); 
    ECP_BN462_copy(&c_q, &Q); 
    ECP_BN462_add(&c_q, &temp1);
    ECP_BN462_copy(&G, &PK.g); 
    ECP2_BN462_copy(&H_hat, &PK.h_hat); 
    ECP2_BN462_copy(&H_theta, &H_hat); 
    PAIR_BN462_G2mul(&H_theta, theta); 

    for (i = 0; i < 2; ++i){
        BIG_BN462_randtrunc(R_M[i], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(R_R[i], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(R_Q[i], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(R_G[i], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(R_H[i], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(R_theta[i], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
    }
    commit_G1( C_M,  R_M, M,  PK.crs_u);
    commit_G1( C_R,  R_R, R,  PK.crs_u);
    commit_G1( C_Q,  R_Q, Q,  PK.crs_u);
    commit_G1( C_G,  R_G, G,  PK.crs_u);
    commit_G2( C_H,  R_H, H_hat,  PK.crs_v);
    commit_G2( C_theta,  R_theta, H_theta,  PK.crs_v);

// 1. 2. PROVE THAT CT IS WELL-FORMED: EQUATIONS 3A, 3B
    ECP_BN462_copy(&g1_inv, &PK.g1);
    ECP_BN462_neg(&g1_inv); 
    ECP_BN462_copy(&h1_inv, &PK.h1);
    ECP_BN462_neg(&h1_inv); 

    ECP_BN462 A5a[2] ={c_1, g1_inv}, pi_eq5a[2], A5b[2] ={c_2, h1_inv}, pi_eq5b[2];
    BIG_BN462 S5a[2][2];

    for (i = 0; i < 2; ++i){
        BIG_464_60_rcopy(S5a[0][i], R_H[i]);
        BIG_464_60_rcopy(S5a[1][i], R_theta[i]);
    }

    lenA = 2;
    prove_pairing_linear_AY(lenA, pi_eq5a, A5a, S5a); 
    prove_pairing_linear_AY(lenA, pi_eq5b, A5b, S5a); 

// 3. PROVE THAT CT IS WELL-FORMED: EQUATION 3C
    ECP_BN462 X[1], qu_temp, qu_temp1, f_inv, A[2], C_f1[2], C_f2[2], C_f3[2], C[1][2], theta_eq5c[2][2];
    ECP2_BN462 Y[2], B[1], D[2][2], pi_eq5c[2][2];
    BIG_BN462 Gamma[1][2], R_f1[2], R_f2[2], R_f3[2], R_5[1][2], T[2][2];

    ECP_BN462_copy(&qu_temp, &PK.f1); 
    PAIR_BN462_G1mul(&qu_temp, theta);
    ECP_BN462_copy(&X[0], &qu_temp); 
    ECP2_BN462_copy(&Y[0], &H_hat); 
    ECP2_BN462_copy(&Y[1], &H_theta);
    ECP_BN462_copy(&f_inv, &PK.f1);
    ECP_BN462_neg(&f_inv);
    iota1(f_inv, A); 
    ECP2_BN462_inf(&B[0]); 
    BIG_464_60_one(Gamma[0][0]); 
    BIG_464_60_zero(Gamma[0][1]);
                                 
    for(i =0; i<2; i++){ 
        BIG_BN462_modneg(R_f1[i], R_M[i], PK.r); 
        BIG_BN462_modneg(R_f2[i], R_R[i], PK.r);
        BIG_BN462_modneg(R_f3[i], R_Q[i], PK.r);
        BIG_464_60_rcopy(R_5[0][i], R_f1[i]); 
        for(j =0; j<2; j++){
            BIG_BN462_randtrunc(T[i][j], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
        }
    }

    lenA = 2, lenB = 1, lenX = 1, lenY = 2; 
    prove_pairing_quadratic(lenA, lenB, lenX, lenY, pi_eq5c, theta_eq5c, PK.crs_u, PK.crs_v, X, Y, A, B, Gamma, R_5, S5a, T); 

//4. PROVE THAT CT IS WELL-FORMED: EQUATIONS 3D
    ECP2_BN462 pi_eq5d[2][2];
    ECP_BN462 theta_eq5d[2][2];

    ECP_BN462_copy(&qu_temp, &PK.f2); 
    PAIR_BN462_G1mul(&qu_temp, theta);
    ECP_BN462_copy(&X[0], &qu_temp);
    ECP_BN462_copy(&f_inv, &PK.f2); 
    ECP_BN462_neg(&f_inv);
    iota1(f_inv, A); 
                                 
    for(i =0; i<2; i++){ 
        BIG_464_60_rcopy(R_5[0][i], R_f2[i]);
        for(j =0; j<2; j++){
            BIG_BN462_randtrunc(T[i][j], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
        }
    }

    lenA = 2, lenB = 1, lenX = 1, lenY = 2; 
    prove_pairing_quadratic(lenA, lenB, lenX, lenY, pi_eq5d, theta_eq5d, PK.crs_u, PK.crs_v, X, Y, A, B, Gamma, R_5, S5a, T); 

//5. PROVE THAT CT IS WELL-FORMED: EQUATION 3E
    ECP2_BN462 pi_eq5e[2][2];
    ECP_BN462 theta_eq5e[2][2];

    ECP_BN462_copy(&qu_temp, &PK.f3); 
    PAIR_BN462_G1mul(&qu_temp, theta);
    ECP_BN462_copy(&X[0], &qu_temp);
    ECP_BN462_copy(&f_inv, &PK.f3); 
    ECP_BN462_neg(&f_inv);
    iota1(f_inv, A); 
                                      
    for(i =0; i<2; i++){ 
        BIG_464_60_rcopy(R_5[0][i], R_f3[i]); 
        for(j =0; j<2; j++){
            BIG_BN462_randtrunc(T[i][j], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
        }
    }
    lenA = 2, lenB = 1, lenX = 1, lenY = 2;
    prove_pairing_quadratic(lenA, lenB, lenX, lenY, pi_eq5e, theta_eq5e, PK.crs_u, PK.crs_v, X, Y, A, B, Gamma, R_5, S5a, T);

//6. PROVE PROVE THAT G AND H_hat ARE WELL-FORMED: EQUATION 2
    ECP2_BN462 Y_eq4[1], pi_eq4[2][2], D_4[1][2];
    ECP_BN462 A_eq4[1], theta_eq4[2][2];
    BIG_BN462 Gamma_eq4[1][2], R_4[1][2], S_4[1][2];

    ECP_BN462_copy(&X[0], &G); 
    ECP2_BN462_copy(&Y_eq4[0], &H_hat); 
    ECP_BN462_copy(&A_eq4[0], &PK.g); 
    ECP2_BN462_copy(&B[0], &PK.h_hat); 
    ECP2_BN462_neg(&B[0]);

    for(i =0; i<2; i++){ 
        BIG_464_60_zero(Gamma_eq4[0][i]);
        BIG_464_60_rcopy(R_4[0][i], R_G[i]);
        BIG_464_60_rcopy(S_4[0][i], R_H[i]); 
        for(j =0; j<2; j++){
            BIG_BN462_randtrunc(T[i][j], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
        }
    }

    lenA = 1, lenB = 1, lenX = 1, lenY = 1;
    prove_pairing_quadratic(lenA, lenB, lenX, lenY, pi_eq4, theta_eq4, PK.crs_u, PK.crs_v, X, Y_eq4, A_eq4, B, Gamma_eq4, R_4, S_4, T); // len A B, X, Y = 2112
    
//7. 8. COMPUTE THE PROOF OF OPENINGS: EQUATIONS 4A 4B
    ECP2_BN462 d1_hat_inv, d2_hat_inv, pi1_open[2], pi2_open[2];
    ECP2_BN462_copy(&d1_hat_inv, &d1_hat);
    ECP2_BN462_neg(&d1_hat_inv);
    ECP2_BN462_copy(&d2_hat_inv, &d2_hat);
    ECP2_BN462_neg(&d2_hat_inv);  

    ECP2_BN462 B1open[4] ={PK.g_hat, PK.g1_hat, PK.h1_hat, d1_hat_inv};
    ECP2_BN462 B2open[4] ={PK.g2_hat, PK.h2_hat, d2_hat_inv};
    BIG_BN462 R1open[4][2], R2open[3][2];

    for(i =0; i<2; i++){
        BIG_464_60_rcopy(R1open[0][i], R_M[i]);
        BIG_464_60_rcopy(R1open[1][i], R_R[i]);
        BIG_464_60_rcopy(R1open[2][i], R_Q[i]);
        BIG_464_60_rcopy(R1open[3][i], R_G[i]);
        BIG_464_60_rcopy(R2open[0][i], R_R[i]);
        BIG_464_60_rcopy(R2open[1][i], R_Q[i]);
        BIG_464_60_rcopy(R2open[2][i], R_G[i]);
    }
    
    lenC = 4;
    prove_pairing_linear_XB(lenC, pi1_open, B1open, R1open);
    lenC = 3;
    prove_pairing_linear_XB(lenC, pi2_open, B2open, R2open);

// 9. COMPUTE THE PROOF OF TRACEABILITY: EQUATION 6
    int slen=3;
    BIG_BN462 chi_trace[slen], gamma_trace[slen], R_Z_hat[2], R_R_hat[2], S_sig[2][2];
    ECP2_BN462  Z1[1], R1[1], Z2[1], R2[1], Z3[1], R3[1], msg1[slen], msg2[slen], msg3[slen], C_Z_hat[2], C_R_hat[2];
    ECP_BN462 pk_sig[slen], temp1_i[slen], A_sig[2] ={PK.g1, PK.h1}, pi_sig[2];

    for (i = 0; i < slen; ++i){
        BIG_BN462_randtrunc(chi_trace[i], PK.r, 2 * CURVE_SECURITY_BN462, RNG); 
        BIG_BN462_randtrunc(gamma_trace[i], PK.r, 2 * CURVE_SECURITY_BN462, RNG); 
        ECP_BN462_copy(&pk_sig[i], &PK.g1); 
        ECP_BN462_copy(&temp1_i[i], &PK.h1); 
        PAIR_BN462_G1mul(&temp1_i[i], gamma_trace[i]); 
        PAIR_BN462_G1mul(&pk_sig[i], chi_trace[i]); 
        ECP_BN462_add(&pk_sig[i], &temp1_i[i]); 
    }
    // sign row 1
    ECP2_BN462_copy(&msg1[0], &PK.g_hat);
    ECP2_BN462_copy(&msg1[1], &d1_hat);
    ECP2_BN462_copy(&msg1[2], &d2_hat);
    ECP2_BN462_inf(&Z1[0]);
    ECP2_BN462_inf(&R1[0]);
    sign_lhspG2(slen, Z1,  R1,  gamma_trace,  chi_trace, msg1); 

    // sign row 2 
    ECP2_BN462_inf(&msg2[0]);
    ECP2_BN462_copy(&msg2[1], &PK.g1_hat);
    ECP2_BN462_copy(&msg2[2], &PK.g2_hat);
    ECP2_BN462_inf(&Z2[0]);
    ECP2_BN462_inf(&R2[0]);
    sign_lhspG2(slen, Z2,  R2,  gamma_trace,  chi_trace, msg2);

    // sign row 3 
    ECP2_BN462_inf(&msg3[0]);
    ECP2_BN462_copy(&msg3[1], &PK.h1_hat);
    ECP2_BN462_copy(&msg3[2], &PK.h2_hat);
    ECP2_BN462_inf(&Z3[0]);
    ECP2_BN462_inf(&R3[0]);
    sign_lhspG2(slen, Z3,  R3,  gamma_trace,  chi_trace, msg3);
    
    for (i = 0; i < 2; ++i){
        BIG_BN462_randtrunc(R_Z_hat[i], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(R_R_hat[i], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
    }
    ECP2_BN462 Z_sig1, R_sig1; 
    ECP2_BN462_copy(&Z_sig1, &Z1[0]);
    ECP2_BN462_copy(&R_sig1, &R1[0]);
    commit_G2( C_Z_hat,  R_Z_hat, Z_sig1,  PK.crs_v_prime); 
    commit_G2( C_R_hat,  R_R_hat, R_sig1,  PK.crs_v_prime);
    
    for (i = 0; i < 2; ++i){
        BIG_464_60_rcopy(S_sig[0][i], R_Z_hat[i]);
        BIG_464_60_rcopy(S_sig[1][i], R_R_hat[i]);
    }
    
    lenA = 2;
    prove_pairing_linear_AY(lenA, pi_sig, A_sig, S_sig); 

//10. COMPUTE THE TCCA PROOF: EQUATION 7
    ECP_BN462 tau, C_X[2], A_inf, B_inf, C_A[2], C_B[2];
    ECP2_BN462 pk_sxdh,  pi_ss[2], B_ss[3] ={PK.g_hat, PK.h_hat, pk_sxdh};
    BIG_BN462 R_A[2], R_B[2], R_X[2], tau_big, R_ss[3][2];
    
    char oct[128];
    octet  tau_oct = {1, sizeof(PK.r), oct}; 

    ECP_BN462_inf(&tau);
    ECP_BN462_inf(&A_inf);
    ECP_BN462_inf(&B_inf);
    derive_commitment1(C_X,  C_G,  PK.g); 

    for (i = 0; i < 2; ++i){
        BIG_BN462_modneg(R_X[i], R_G[i], PK.r);
        BIG_BN462_randtrunc(R_A[i], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
        BIG_BN462_randtrunc(R_B[i], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
    }
    commit_G1( C_A,  R_A, A_inf,  PK.crs_u);
    commit_G1( C_B,  R_B, B_inf,  PK.crs_u);
    
    for (i = 0; i < 3; ++i){
        ECP_BN462_add(&tau, &pk_sig[i]);
    }
    
    ECP_BN462_toOctet(&tau_oct, &tau, true);
    BIG_464_60_fromBytes(tau_big, tau_oct.val);
    ECP2_BN462_copy(&pk_sxdh, &PK.f_hat[1]);
    PAIR_BN462_G2mul(&pk_sxdh, tau_big);
    ECP2_BN462_add(&pk_sxdh, &PK.f_hat[0]);
    ECP2_BN462_neg(&pk_sxdh);
    ECP2_BN462_copy(&B_ss[2], &pk_sxdh);    
    
    for (i = 0; i < 2; ++i){
        BIG_464_60_rcopy(R_ss[0][i], R_A[i]);
        BIG_464_60_rcopy(R_ss[1][i], R_B[i]);
        BIG_464_60_rcopy(R_ss[2][i], R_X[i]);
    }
    
    lenC = 3;
    prove_pairing_linear_XB(lenC, pi_ss, B_ss, R_ss);

//11. COMPUTE THE WELL-FORMEDNESS PROOF OF THE VOTE: EQUATION 8A
    ECP2_BN462  C_M_hat[2], M_hat, Y_10a[1], pi_10a[2][2], D_10a[1][2];
    ECP_BN462 A_10a[1], theta_10a[2][2];
    BIG_BN462 Gamma_10a[1][2], R_10a[1][2], S_10a[1][2], R_M_hat[2];

    ECP2_BN462_copy(&M_hat, &PK.g_hat); 
    PAIR_BN462_G2mul(&M_hat, m);
    for (i = 0; i < 2; ++i){
        BIG_BN462_randtrunc(R_M_hat[i], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
    }
    commit_G2( C_M_hat,  R_M_hat, M_hat,  PK.crs_v);

    ECP_BN462_copy(&X[0], &M); 
    ECP2_BN462_copy(&Y_10a[0], &M_hat); 
    ECP_BN462_copy(&A_10a[0], &PK.g); 
    ECP2_BN462_copy(&B[0], &PK.g_hat); 
    ECP2_BN462_neg(&B[0]);
    
    for(i =0; i<2; i++){ 
        BIG_464_60_zero(Gamma_10a[0][i]);
        BIG_464_60_rcopy(R_10a[0][i], R_M[i]);
        BIG_464_60_rcopy(S_10a[0][i], R_M_hat[i]); 
        for(j =0; j<2; j++){
            BIG_BN462_randtrunc(T[i][j], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
        }
    }

    lenA = 1, lenB = 1, lenX = 1, lenY = 1;
    prove_pairing_quadratic(lenA, lenB, lenX, lenY, pi_10a, theta_10a, PK.crs_u, PK.crs_v, X, Y_10a, A_10a, B, Gamma_10a, R_10a, S_10a, T);
    
// 12. COMPUTE THE WELL-FORMEDNESS PROOF OF THE VOTE: EQUATION 8B
    ECP_BN462  A_10b[1], theta_10b[2][2];
    ECP2_BN462 Y_10b[1], pi_10b[2][2];
    BIG_BN462 Gamma_10b[1][1], R_M_hat_inv[2], R_10b[1][2], S_10b[1][2];

    ECP_BN462_copy(&X[0], &M); 
    ECP2_BN462_copy(&Y_10b[0], &PK.g_hat); 
    PAIR_BN462_G2mul(&Y_10b[0], m);
    ECP2_BN462_neg(&Y_10b[0]);
    ECP2_BN462_add(&Y_10b[0], &PK.g_hat); 
    ECP2_BN462_inf(&B[0]); 
    ECP_BN462_inf(&A_10b[0]); 
    BIG_464_60_one(Gamma_10b[0][0]); 
                                    
    for(i =0; i<2; i++){ 
        BIG_BN462_modneg(R_M_hat_inv[i], R_M_hat[i], PK.r); 
        BIG_464_60_rcopy(R_10b[0][i], R_M[i]);
        BIG_464_60_rcopy(S_10b[0][i], R_M_hat_inv[i]); 
        for(j =0; j<2; j++){
            BIG_BN462_randtrunc(T[i][j], PK.r, 2 * CURVE_SECURITY_BN462, RNG);
        }
    }

    lenA = 1, lenB = 1, lenX = 1, lenY = 1; 
    prove_pairing_quadratic(lenA, lenB, lenX, lenY, pi_10b, theta_10b, PK.crs_u, PK.crs_v, X, Y_10b, A_10b, B, Gamma_10b, R_10b, S_10b, T); 


    // Set values to ballot
    for(i =0; i<2; i++){
        ECP_BN462_copy(&ballot->C_M[i], &C_M[i]);
        ECP_BN462_copy(&ballot->C_R[i], &C_R[i]);
        ECP_BN462_copy(&ballot->C_Q[i], &C_Q[i]);
        ECP_BN462_copy(&ballot->C_G[i], &C_G[i]);
        ECP_BN462_copy(&ballot->C_A[i], &C_A[i]);
        ECP_BN462_copy(&ballot->C_B[i], &C_B[i]);
        ECP2_BN462_copy(&ballot->C_H[i], &C_H[i]);
        ECP2_BN462_copy(&ballot->C_theta[i], &C_theta[i]);
        ECP2_BN462_copy(&ballot->C_Z_hat[i], &C_Z_hat[i]);
        ECP2_BN462_copy(&ballot->C_R_hat[i], &C_R_hat[i]);
        ECP2_BN462_copy(&ballot->C_M_hat[i], &C_M_hat[i]);
        ECP_BN462_copy(&ballot->pi_eq5a[i], &pi_eq5a[i]);
        ECP_BN462_copy(&ballot->pi_eq5b[i], &pi_eq5b[i]);
        ECP2_BN462_copy(&ballot->pi1_open[i], &pi1_open[i]);
        ECP2_BN462_copy(&ballot->pi2_open[i], &pi2_open[i]);
        ECP_BN462_copy(&ballot->pi_sig[i], &pi_sig[i]);
        ECP2_BN462_copy(&ballot->pi_ss[i], &pi_ss[i]);
        

        for(j =0; j<2; j++){
            ECP2_BN462_copy(&ballot->pi_eq5c[i][j], &pi_eq5c[i][j]);
            ECP_BN462_copy(&ballot->theta_eq5c[i][j], &theta_eq5c[i][j]);
            ECP2_BN462_copy(&ballot->pi_eq5d[i][j], &pi_eq5d[i][j]);
            ECP_BN462_copy(&ballot->theta_eq5d[i][j], &theta_eq5d[i][j]);
            ECP2_BN462_copy(&ballot->pi_eq5e[i][j], &pi_eq5e[i][j]);
            ECP_BN462_copy(&ballot->theta_eq5e[i][j], &theta_eq5e[i][j]);
            ECP2_BN462_copy(&ballot->pi_eq4[i][j], &pi_eq4[i][j]);
            ECP_BN462_copy(&ballot->theta_eq4[i][j], &theta_eq4[i][j]);
            ECP2_BN462_copy(&ballot->pi_10a[i][j], &pi_10a[i][j]);
            ECP_BN462_copy(&ballot->theta_10a[i][j], &theta_10a[i][j]);
            ECP2_BN462_copy(&ballot->pi_10b[i][j], &pi_10b[i][j]);
            ECP_BN462_copy(&ballot->theta_10b[i][j], &theta_10b[i][j]);
        }
        
    }
    for(i =0; i<3; i++){
        ECP_BN462_copy(&ballot->pk_sig[i], &pk_sig[i]);
    }
    ECP2_BN462_copy(&ballot->Z2[0], &Z2[0]);
    ECP2_BN462_copy(&ballot->R2[0], &R2[0]);
    ECP2_BN462_copy(&ballot->Z3[0], &Z3[0]);
    ECP2_BN462_copy(&ballot->R3[0], &R3[0]);
    ECP_BN462_copy(&ballot->c_1, &c_1); 
    ECP_BN462_copy(&ballot->c_2, &c_2); 
    ECP_BN462_copy(&ballot->c_0, &c_0); 
    ECP_BN462_copy(&ballot->c_r, &c_r); 
    ECP_BN462_copy(&ballot->c_q, &c_q); 
    ECP2_BN462_copy(&ballot->d1_hat, &d1_hat); 
    ECP2_BN462_copy(&ballot->d2_hat, &d2_hat); 
}

bool verify(pk_simple_ballot PK, simple_ballot ballot){ 
    bool is_valid[12]; 
    int i = 0, j = 0, lenA = 0, lenB = 0, lenC = 0, lenD = 0, lenX = 0, lenY = 0;

//1. 2. VERIFY EQUATIONS 3A, 3B
    ECP2_BN462 D5a[2][2]; 
    ECP_BN462 g1_inv, h1_inv;
    FP12_BN462 t_eq5a, t_eq5b, t1_open, t2_open;
    lenA = 2;

    ECP_BN462_copy(&g1_inv, &PK.g1);
    ECP_BN462_neg(&g1_inv); 
    ECP_BN462_copy(&h1_inv, &PK.h1);
    ECP_BN462_neg(&h1_inv); 
    ECP_BN462 A5a[2] ={ballot.c_1, g1_inv}, A5b[2] ={ballot.c_2, h1_inv};

    for(i =0; i<2; i++){
        ECP2_BN462_copy(&D5a[0][i], &ballot.C_H[i]);
        ECP2_BN462_copy(&D5a[1][i], &ballot.C_theta[i]);
    }
   
    FP12_BN462_one(&t_eq5a);
    is_valid[0] = verify_pairing_linear_AY(lenA, PK.crs_v, A5a, D5a, ballot.pi_eq5a, t_eq5a);
    FP12_BN462_one(&t_eq5b);
    is_valid[1] = verify_pairing_linear_AY(lenA, PK.crs_v, A5b, D5a, ballot.pi_eq5b, t_eq5b);

//3. VERIFY EQUATION 3C
    FP12_BN462 t_eq5c;
    ECP_BN462  C_f1[2], C_f2[2], C_f3[2], C[1][2],  A[2], f_inv;
    ECP2_BN462 D[2][2], B[1];
    BIG_BN462 Gamma[1][2];

    ECP_BN462_copy(&f_inv, &PK.f1);
    ECP_BN462_neg(&f_inv);
    iota1(f_inv, A); 
    ECP2_BN462_inf(&B[0]); 
    BIG_464_60_one(Gamma[0][0]); 
    BIG_464_60_zero(Gamma[0][1]);
    FP12_BN462_one(&t_eq5c);
    derive_commitment1(C_f1,  ballot.C_M,  ballot.c_0);

    for(i =0; i<2; i++){
        ECP_BN462_copy(&C[0][i], &C_f1[i]);
        ECP2_BN462_copy(&D[0][i], &ballot.C_H[i]);
        ECP2_BN462_copy(&D[1][i], &ballot.C_theta[i]);
    }
    
    lenA =2, lenC = 1, lenD = 2;
    is_valid[2]=  verify_pairing_quadratic1(lenA, lenC, lenD, PK.crs_u, PK.crs_v, C, D, A, B, Gamma, ballot.pi_eq5c, ballot.theta_eq5c, t_eq5c); 
    
//4. VERIFY EQUATION 3D
    FP12_BN462 t_eq5d;

    ECP_BN462_copy(&f_inv, &PK.f2); 
    ECP_BN462_neg(&f_inv);
    iota1(f_inv, A); 
    FP12_BN462_one(&t_eq5d);
    derive_commitment1(C_f2,  ballot.C_R,  ballot.c_r); 

    for(i =0; i<2; i++){
        ECP_BN462_copy(&C[0][i], &C_f2[i]);
    }
    lenA =2, lenC = 1, lenD = 2;
    is_valid[3] = verify_pairing_quadratic1(lenA, lenC, lenD, PK.crs_u, PK.crs_v, C, D, A, B, Gamma, ballot.pi_eq5d, ballot.theta_eq5d, t_eq5d);
    
//5. VERIFY EQUATION 3E
    FP12_BN462 t_eq5e;

    ECP_BN462_copy(&f_inv, &PK.f3); 
    ECP_BN462_neg(&f_inv);
    iota1(f_inv, A); 
    FP12_BN462_one(&t_eq5e);
    derive_commitment1(C_f3,  ballot.C_Q,  ballot.c_q); 

    for(i =0; i<2; i++){
        ECP_BN462_copy(&C[0][i], &C_f3[i]);
    }
    lenA =2, lenC = 1, lenD = 2;
    is_valid[4] = verify_pairing_quadratic1(lenA, lenC, lenD, PK.crs_u, PK.crs_v, C, D, A, B, Gamma, ballot.pi_eq5e, ballot.theta_eq5e, t_eq5e); 
    
// 6. VERIFY EQUATION 2
    FP12_BN462 t_eq4;
    ECP2_BN462 Y_eq4[1], D_4[1][2];
    ECP_BN462 A_eq4[1];
    BIG_BN462 Gamma_eq4[1][2];

    FP12_BN462_one(&t_eq4);
    ECP_BN462_copy(&A_eq4[0], &PK.g); 
    ECP2_BN462_copy(&B[0], &PK.h_hat); 
    ECP2_BN462_neg(&B[0]);

    for(i =0; i<2; i++){
        BIG_464_60_zero(Gamma_eq4[0][i]); 
        ECP_BN462_copy(&C[0][i], &ballot.C_G[i]);
        ECP2_BN462_copy(&D_4[0][i], &ballot.C_H[i]);
    }
    lenA = 1, lenC = 1, lenD = 1;
    is_valid[5] = verify_pairing_quadratic(lenA, lenC, lenD, PK.crs_u, PK.crs_v, C, D_4, A_eq4, B, Gamma_eq4, ballot.pi_eq4, ballot.theta_eq4, t_eq4);

//7. 8. VERIFY EQUATIONS 4A 4B
    ECP_BN462 Copen1[4][2], Copen2[3][2]; 
    ECP2_BN462 d1_hat_inv, d2_hat_inv;

    ECP2_BN462_copy(&d1_hat_inv, &ballot.d1_hat);
    ECP2_BN462_neg(&d1_hat_inv);
    ECP2_BN462_copy(&d2_hat_inv, &ballot.d2_hat);
    ECP2_BN462_neg(&d2_hat_inv); 
    ECP2_BN462 B1open[4] ={PK.g_hat, PK.g1_hat, PK.h1_hat, d1_hat_inv};
    ECP2_BN462 B2open[4] ={PK.g2_hat, PK.h2_hat, d2_hat_inv};

    for(i =0; i<2; i++){
        ECP_BN462_copy(&Copen1[0][i], &ballot.C_M[i]);
        ECP_BN462_copy(&Copen1[1][i], &ballot.C_R[i]);
        ECP_BN462_copy(&Copen1[2][i], &ballot.C_Q[i]);
        ECP_BN462_copy(&Copen1[3][i], &ballot.C_G[i]);

        ECP_BN462_copy(&Copen2[0][i], &ballot.C_R[i]);
        ECP_BN462_copy(&Copen2[1][i], &ballot.C_Q[i]);
        ECP_BN462_copy(&Copen2[2][i], &ballot.C_G[i]);
    }

    lenC = 4;
    FP12_BN462_one(&t1_open);
    is_valid[6] = verify_pairing_linear_XB(lenC, PK.crs_u, B1open, Copen1, ballot.pi1_open, t1_open);

    lenC = 3;
    FP12_BN462_one(&t2_open);
    is_valid[7] = verify_pairing_linear_XB(lenC, PK.crs_u, B2open, Copen2, ballot.pi2_open, t2_open);
    
//9. VERIFY EQUATION 6 
    int slen = 3;
    FP12_BN462 t_sig, temp_sig;
    ECP2_BN462 D_sig[2][2], msg1[slen]; 
    ECP_BN462  A_sig[2] ={PK.g1, PK.h1};
    ECP2_BN462_copy(&msg1[0], &PK.g_hat);
    ECP2_BN462_copy(&msg1[1], &ballot.d1_hat);
    ECP2_BN462_copy(&msg1[2], &ballot.d2_hat);

    for(i =0; i<2; i++){
        ECP2_BN462_copy(&D_sig[0][i], &ballot.C_Z_hat[i]);
        ECP2_BN462_copy(&D_sig[1][i], &ballot.C_R_hat[i]);
    }
    
    FP12_BN462_one(&t_sig);
    for (i = 0; i < slen; ++i){
        temp_sig= PAIR_BN462_e(msg1[i], ballot.pk_sig[i]);
        FP12_BN462_ssmul(&t_sig,&temp_sig);
    }
    lenA = 2;
    is_valid[8] = verify_pairing_linear_AY(lenA, PK.crs_v_prime, A_sig, D_sig, ballot.pi_sig, t_sig);

//10. VERIFY EQUATION 7
    FP12_BN462 t_ss;
    ECP_BN462 C_ss[3][2], C_X[2], tau; 
    ECP2_BN462 pk_sxdh;
    BIG_BN462 tau_big;
    char oct[128];
    octet  tau_oct = {1, sizeof(PK.r), oct};

    ECP_BN462_inf(&tau);
    derive_commitment1(C_X,  ballot.C_G,  PK.g); 
    for(i =0; i<2; i++){
        ECP_BN462_copy(&C_ss[0][i], &ballot.C_A[i]);
        ECP_BN462_copy(&C_ss[1][i], &ballot.C_B[i]);
        ECP_BN462_copy(&C_ss[2][i], &C_X[i]);
    }

    for (i = 0; i < 3; ++i){
        ECP_BN462_add(&tau, &ballot.pk_sig[i]);
    }
    ECP_BN462_toOctet(&tau_oct, &tau, true);
    BIG_464_60_fromBytes(tau_big, tau_oct.val);
    ECP2_BN462_copy(&pk_sxdh, &PK.f_hat[1]);
    PAIR_BN462_G2mul(&pk_sxdh, tau_big);
    ECP2_BN462_add(&pk_sxdh, &PK.f_hat[0]);
    ECP2_BN462_neg(&pk_sxdh);
    ECP2_BN462 B_ss[4] ={PK.g_hat, PK.h_hat, pk_sxdh};
    FP12_BN462_one(&t_ss);

    lenC = 3;
    is_valid[9] = verify_pairing_linear_XB(lenC, PK.crs_u, B_ss, C_ss, ballot.pi_ss, t_ss);

//11. VERIFY EQUATION 8A
    FP12_BN462 t_10a;
    ECP2_BN462 Y_10a[1], D_10a[1][2];
    ECP_BN462 A_10a[1];
    BIG_BN462 Gamma_10a[1][2];

    FP12_BN462_one(&t_10a);
    ECP_BN462_copy(&A_10a[0], &PK.g); 
    ECP2_BN462_copy(&B[0], &PK.g_hat);
    ECP2_BN462_neg(&B[0]);

    for (i = 0; i < 2; ++i){
        BIG_464_60_zero(Gamma_10a[0][i]); 
        ECP_BN462_copy(&C[0][i], &ballot.C_M[i]);
        ECP2_BN462_copy(&D_10a[0][i], &ballot.C_M_hat[i]);
    }
    
    lenA = 1, lenC = 1, lenD = 1;
    is_valid[10] = verify_pairing_quadratic(lenA, lenC, lenD, PK.crs_u, PK.crs_v, C, D_10a, A_10a, B, Gamma_10a, ballot.pi_10a, ballot.theta_10a, t_10a); 
    
//12. VERIFY EQUATION 8B
    FP12_BN462 t_10b;
    ECP2_BN462  C_M_hat_inv[2], D_10b[1][2];
    ECP_BN462  A_10b[1];
    BIG_BN462 Gamma_10b[1][1];
    
    ECP2_BN462_inf(&B[0]); 
    ECP_BN462_inf(&A_10b[0]); 
    BIG_464_60_one(Gamma_10b[0][0]); 
    FP12_BN462_one(&t_10b);
    derive_commitment2(C_M_hat_inv,  ballot.C_M_hat,  PK.g_hat);

    for (i = 0; i < 2; ++i){
        ECP_BN462_copy(&C[0][i], &ballot.C_M[i]);
        ECP2_BN462_copy(&D_10b[0][i], &C_M_hat_inv[i]);
    }
    lenA = 1, lenC = 1, lenD = 1;
    is_valid[11] = verify_pairing_quadratic1(lenA, lenC, lenD, PK.crs_u, PK.crs_v, C, D_10b, A_10b, B, Gamma_10b, ballot.pi_10b, ballot.theta_10b, t_10b);
     
    for ( i = 0; i < 12; i++){
        if(is_valid[i] == false)
            return false;
    } 
    return true;
}

void decrypt(BIG_BN462 decrypted_m, sk_simple_ballot SK, pk_simple_ballot PK, simple_ballot ballot){
    BIG_BN462 x, y, unit;
    BIG_464_60_zero(x);
    BIG_464_60_one(unit);
    
    if(! verify(PK, ballot)){
        printf("Verification returns 0!\n");
    } 
    else{
        int count = 0;
        ECP_BN462 M, temp;
        ECP_BN462_copy(&M,&ballot.c_0);
        ECP_BN462_copy(&temp,&ballot.c_1);
        PAIR_BN462_G1mul(&temp,SK.alpha[0]);
        ECP_BN462_neg(&temp);
        ECP_BN462_add(&M,&temp);
        ECP_BN462_copy(&temp,&ballot.c_2);
        PAIR_BN462_G1mul(&temp,SK.beta[0]);
        ECP_BN462_neg(&temp);
        ECP_BN462_add(&M,&temp);
        BIG_464_60_modneg(x, unit, PK.r); 
        do {
            BIG_464_60_modadd(x, x, unit, PK.r); 
            ECP_BN462_copy(&temp,&PK.g);
            PAIR_BN462_G1mul(&temp, x); 
        } 
        while(ECP_BN462_equals(&temp,&M) != 1); 
        BIG_464_60_rcopy(decrypted_m, x);
    }
    
}
