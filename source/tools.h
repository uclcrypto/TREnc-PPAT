#include "../library/pair_BN462.h"

#define BIG_BN462 BIG_464_60
#define BIG_BN462_rcopy BIG_464_60_rcopy
#define BIG_BN462_randomnum BIG_464_60_randomnum
#define BIG_BN462_randtrunc BIG_464_60_randtrunc
#define BIG_BN462_modmul BIG_464_60_modmul
#define BIG_BN462_modadd BIG_464_60_modadd
#define BIG_BN462_modneg BIG_464_60_modneg
#define BIG_BN462_output BIG_464_60_output

// complex ballots
typedef struct pk_complex_ballot {
    ECP_BN462 g, h, g1, h1, f, k_i[2], crs_u[2][2], crs_u_prime[2][2];
    ECP2_BN462  g_hat, h_hat, f1_hat, f2_hat, crs_v[2][2];
    BIG_BN462 order;
} pk_complex_ballot;

typedef struct sk_complex_ballot {
    BIG_BN462 a[3], b[3];
} sk_complex_ballot;

typedef struct complex_ballot {
    ECP_BN462 C_G[2], C_M[2], C_G_gamma[2], C_Z1[2], C_R1[2], cm_0, cm_1, cm_2, d1, d2;
    ECP2_BN462 crhat_0, crhat_1, crhat_2, cshat, pk_sig[3], C_G_hat[2], C_theta[2], C_R_hat[2], C_S_hat[2], C_A[2], C_B[2];
    ECP_BN462 pi_eq5a[2], pi_eq5b[2], pi_12b[2], pi_ss[2];
    ECP2_BN462 pi_eq5e[2][2], pi_11c[2], pi_11d[2], pi_eq4[2][2], pi_12a[2][2], pi_14[2], pi_11f[2][2], pi_11g[2][2];
    ECP_BN462 theta_eq5e[2][2], theta_eq4[2][2], theta_12a[2][2], theta_11f[2][2], theta_11g[2][2];
}complex_ballot;


// simple ballots
typedef struct pk_simple_ballot {
    ECP_BN462 g, g1, h1, f1, f2, f3, crs_u[2][2];;
    ECP2_BN462  g_hat, h_hat, g1_hat, g2_hat, h1_hat, h2_hat, crs_v[2][2], f_hat[2], crs_v_prime[2][2];
    BIG_BN462 r;
} pk_simple_ballot;

typedef struct sk_simple_ballot {
    BIG_BN462 alpha[3], beta[3];
} sk_simple_ballot;

typedef struct simple_ballot {
    ECP_BN462 c_1, c_2, c_0, c_r, c_q, C_M[2], C_R[2], C_Q[2], C_G[2], C_A[2], C_B[2], pk_sig[3];  
    ECP2_BN462 C_H[2], C_theta[2], C_Z_hat[2], C_R_hat[2], Z2[1], R2[1], Z3[1], R3[1], C_M_hat[2];
    ECP2_BN462  d1_hat, d2_hat, pi_eq5c[2][2], pi_eq5d[2][2], pi_eq5e[2][2], pi_eq4[2][2], pi1_open[2], pi2_open[2], pi_ss[2], pi_10a[2][2], pi_10b[2][2];
    ECP_BN462 pi_eq5a[2], pi_eq5b[2], theta_eq5c[2][2], theta_eq5d[2][2], theta_eq5e[2][2], theta_eq4[2][2], pi_sig[2], theta_10a[2][2], theta_10b[2][2];

} simple_ballot;

void iota1(ECP_BN462 x1, ECP_BN462 output[]){
    ECP_BN462_inf(&output[0]); 
    ECP_BN462_copy(&output[1], &x1); 
}

void iota2(ECP2_BN462 x1, ECP2_BN462 output[]){
    ECP2_BN462_inf(&output[0]); 
    ECP2_BN462_copy(&output[1], &x1); 
}

void iotat(FP12_BN462 x1, FP12_BN462 output[2][2]){
    FP12_BN462_one(&output[0][0]); 
    FP12_BN462_one(&output[0][1]);
    FP12_BN462_one(&output[1][0]);
    FP12_BN462_copy(&output[1][1], &x1); 
}

FP12_BN462 PAIR_BN462_e(ECP2_BN462 g2, ECP_BN462 g1){
    FP12_BN462 w;
    PAIR_BN462_ate(&w, &g2, &g1);
    PAIR_BN462_fexp(&w);
    return w;
}

void commit_G1(ECP_BN462 C[], BIG_BN462 R[],ECP_BN462 M, ECP_BN462 crs_u[2][2]){
    ECP_BN462 temp;

    ECP_BN462_copy(&C[0], &crs_u[0][0]);  
    PAIR_BN462_G1mul(&C[0], R[0]); 
    ECP_BN462_copy(&temp, &crs_u[1][0]);  
    PAIR_BN462_G1mul(&temp, R[1]); 
    ECP_BN462_add(&C[0], &temp); 
    ECP_BN462_copy(&C[1], &M); 
   
    ECP_BN462_copy(&temp, &crs_u[0][1]);
    PAIR_BN462_G1mul(&temp, R[0]); 
    ECP_BN462_add(&C[1], &temp); 
    ECP_BN462_copy(&temp, &crs_u[1][1]);
    PAIR_BN462_G1mul(&temp, R[1]);
    ECP_BN462_add(&C[1], &temp); 
}

void commit_G2(ECP2_BN462 C[], BIG_BN462 R[], ECP2_BN462 H, ECP2_BN462 crs_v[2][2]){
    ECP2_BN462 temp;

    ECP2_BN462_copy(&C[0], &crs_v[0][0]);  
    PAIR_BN462_G2mul(&C[0], R[0]); 
    ECP2_BN462_copy(&temp, &crs_v[1][0]);  
    PAIR_BN462_G2mul(&temp, R[1]); 
    ECP2_BN462_add(&C[0], &temp); 
    ECP2_BN462_copy(&C[1], &H); 
   
    ECP2_BN462_copy(&temp, &crs_v[0][1]);
    PAIR_BN462_G2mul(&temp, R[0]); 
    ECP2_BN462_add(&C[1], &temp); 
    ECP2_BN462_copy(&temp, &crs_v[1][1]);
    PAIR_BN462_G2mul(&temp, R[1]); 
    ECP2_BN462_add(&C[1], &temp); 
}

void derive_commitment1(ECP_BN462 C_f[], ECP_BN462 C[], ECP_BN462 c){
    ECP_BN462 out[2];
    iota1(c, out);
    for (int i =0; i<2; i++){
        ECP_BN462_copy(&C_f[i], &C[i]); 
        ECP_BN462_neg(&C_f[i]);
        ECP_BN462_add(&C_f[i], &out[i]);
    }
}

void derive_commitment2(ECP2_BN462 C_f[], ECP2_BN462 C_hat[], ECP2_BN462 c_hat){
    ECP2_BN462 out[2];
    iota2(c_hat, out);
    for (int i =0; i<2; i++){
        ECP2_BN462_copy(&C_f[i], &C_hat[i]); 
        ECP2_BN462_neg(&C_f[i]);
        ECP2_BN462_add(&C_f[i], &out[i]);
    }
}

void sign_lhspG1(int slen, ECP_BN462 Z[], ECP_BN462 R[], BIG_BN462 gamma_i[], BIG_BN462 chi_i[], ECP_BN462 m[]){
    //LHSP signature computation for a message vector in G1
    ECP_BN462 temp;
    for (int i =0; i<slen; i++){
        ECP_BN462_copy(&temp,&m[i]);
        PAIR_BN462_G1mul(&temp,chi_i[i]);
        ECP_BN462_add(&Z[0],&temp); 

        ECP_BN462_copy(&temp,&m[i]);
        PAIR_BN462_G1mul(&temp,gamma_i[i]);
        ECP_BN462_add(&R[0],&temp);
    }
}

void sign_lhspG2(int slen, ECP2_BN462 Z[], ECP2_BN462 R[], BIG_BN462 gamma_i[], BIG_BN462 chi_i[], ECP2_BN462 m[]){
    //LHSP signature computation for a message vector in G2 
    ECP2_BN462 temp;
    for (int i =0; i<slen; i++){
        ECP2_BN462_copy(&temp,&m[i]);
        PAIR_BN462_G2mul(&temp,chi_i[i]);
        ECP2_BN462_add(&Z[0],&temp);

        ECP2_BN462_copy(&temp,&m[i]);
        PAIR_BN462_G2mul(&temp,gamma_i[i]);
        ECP2_BN462_add(&R[0],&temp);
    }
}

void prove_pairing_linear_XB(int lenB, ECP2_BN462 proof[], ECP2_BN462 B[], BIG_BN462 R[4][2]){ 
    // Proof of statements e(x_1, b_1)e(x_2, b_2)...e(x_n, b_n) = t with
    // B = [b_1, ..., b_n] and R = [R_1, ..., R_2] are the randomnesses 
    // used in the commitments of [x_1, ..., x_n]
    int i;
    ECP2_BN462 temp;
    for(i=0; i<lenB ; i++){
        ECP2_BN462_inf(&proof[i]); 
    }
    for(i=0; i<lenB ; i++){
        ECP2_BN462_copy(&temp,&B[i]); 
        PAIR_BN462_G2mul(&temp, R[i][0]); 
        ECP2_BN462_add(&proof[0], &temp);

        ECP2_BN462_copy(&temp,&B[i]); 
        PAIR_BN462_G2mul(&temp, R[i][1]); 
        ECP2_BN462_add(&proof[1], &temp); 
    }
}

void prove_pairing_linear_AY(int lenA, ECP_BN462 proof[], ECP_BN462 A[], BIG_BN462 S[2][2]){
    // Proof of statements e(a_1, y_1)e(a_2, y_2)...e(a_n, y_n) = t with
    // A = [a_1, ..., a_n] and S = [S_1, ..., S_2] are the randomnesses 
    // used in the commitments of [y_1, ..., y_n]
    int i;
    ECP_BN462 temp;
    for(i=0; i<lenA ; i++){
        ECP_BN462_inf(&proof[i]) ;
    }
    for(i=0; i<lenA ; i++){
        ECP_BN462_copy(&temp,&A[i]); 
        PAIR_BN462_G1mul(&temp, S[i][0]); 
        ECP_BN462_add(&proof[0], &temp); 
       

        ECP_BN462_copy(&temp,&A[i]); 
        PAIR_BN462_G1mul(&temp, S[i][1]); 
        ECP_BN462_add(&proof[1], &temp); 
    }
}

void prove_pairing_quadratic(int lenA, int lenB, int lenX, int lenY, ECP2_BN462 proof_pi[2][2], ECP_BN462 proof_theta[2][2], ECP_BN462 crs_u[2][2], ECP2_BN462 crs_v[2][2], ECP_BN462 X[], ECP2_BN462 Y[], ECP_BN462 A[], ECP2_BN462 B[], BIG_BN462 Gamma[lenX][lenY], BIG_BN462 R[lenX][2], BIG_BN462 S[lenY][2], BIG_BN462 T[2][2]){
   int i, j, k, l;
   ECP2_BN462 temp2, out_iota2[2];
   ECP_BN462 temp1, out_iota1[2];
   BIG_BN462 temp, r, RGS[2][2];
   BIG_BN462_rcopy(r, CURVE_Order_BN462);

   for(i =0; i<2; i++){ 
        for(j =0; j<2; j++){
            BIG_BN462_modneg(RGS[i][j], T[j][i], r);
            ECP2_BN462_inf(&proof_pi[i][j]);
            ECP_BN462_inf(&proof_theta[i][j]); 
        }
    }

   for (i =0; i<2; i++){ 
        for (j =0; j<2; j++){ 
            for (k =0; k<lenB; k++){
                iota2(B[k], out_iota2); 
                ECP2_BN462_copy(&temp2,&out_iota2[j]); 
                PAIR_BN462_G2mul(&temp2, R[k][i]);  
                ECP2_BN462_add(&proof_pi[i][j], &temp2); 
            }

            for (k =0; k<lenA; k++){
                iota1(A[k], out_iota1); 
                ECP_BN462_copy(&temp1,&out_iota1[j]); 
                PAIR_BN462_G1mul(&temp1, S[k][i]);  
                ECP_BN462_add(&proof_theta[i][j], &temp1); 
            }
        }
   }

    for (i =0; i<2; i++){ 
        for (j =0; j<2; j++){ 
            for (k =0; k<lenX; k++){
                for (l =0; l<lenY; l++){
                    iota2(Y[l], out_iota2); 
                    ECP2_BN462_copy(&temp2,&out_iota2[j]);
                    BIG_BN462_modmul(temp, Gamma[k][l], R[k][i], r) ;
                    PAIR_BN462_G2mul(&temp2, temp);  
                    ECP2_BN462_add(&proof_pi[i][j], &temp2); 

                    iota1(X[k], out_iota1); 
                    ECP_BN462_copy(&temp1,&out_iota1[j]);
                    BIG_BN462_modmul(temp, Gamma[k][l], S[l][i], r) ;
                    PAIR_BN462_G1mul(&temp1, temp);  
                    ECP_BN462_add(&proof_theta[i][j], &temp1); 

                    BIG_BN462_modmul(temp, Gamma[k][l], S[l][j], r) ;
                    BIG_BN462_modmul(temp, temp, R[k][i], r) ; 
                    BIG_BN462_modadd(RGS[i][j], RGS[i][j], temp, r); 
                }
            }
        }
    }

    for (i =0; i<2; i++){ 
        for (j =0; j<2; j++){ 
            for (k =0; k<2; k++){
                ECP2_BN462_copy(&temp2,&crs_v[k][j]); 
                PAIR_BN462_G2mul(&temp2, RGS[i][k]); 
                ECP2_BN462_add(&proof_pi[i][j], &temp2);

                ECP_BN462_copy(&temp1,&crs_u[k][j]); 
                PAIR_BN462_G1mul(&temp1, T[i][k]); 
                ECP_BN462_add(&proof_theta[i][j], &temp1);
            }
        }
    }
}

void prove_pairing_quadratic2(int lenA, int lenB, int lenX, int lenY, ECP2_BN462 proof_pi[2][2], ECP_BN462 proof_theta[2][2], ECP_BN462 crs_u[2][2], ECP2_BN462 crs_v[2][2], ECP_BN462 X[], ECP2_BN462 Y[], ECP_BN462 A[], ECP2_BN462 B[], BIG_BN462 Gamma[lenX][lenY], BIG_BN462 R[lenX][2], BIG_BN462 S[lenY][2], BIG_BN462 T[2][2]){
   int i, j, k, l;
   ECP2_BN462 temp2, out_iota2[2];
   ECP_BN462 temp1, out_iota1[2];
   BIG_BN462 temp, r, RGS[2][2];
   BIG_BN462_rcopy(r, CURVE_Order_BN462);

   for(i =0; i<2; i++){ 
        for(j =0; j<2; j++){
            BIG_BN462_modneg(RGS[i][j], T[j][i], r);
            ECP2_BN462_inf(&proof_pi[i][j]);
            ECP_BN462_inf(&proof_theta[i][j]); 
        }
    }

   for (i =0; i<2; i++){ 
        for (j =0; j<2; j++){ 
            for (k =0; k<lenB; k++){ 
                iota2(B[k], out_iota2); 
                ECP2_BN462_copy(&temp2,&out_iota2[j]); 
                PAIR_BN462_G2mul(&temp2, R[k][i]);  
                ECP2_BN462_add(&proof_pi[i][j], &temp2); 
            }

            for (k =0; k<lenA; k++){
                iota1(A[k], out_iota1); 
                ECP_BN462_copy(&temp1,&out_iota1[j]); 
                PAIR_BN462_G1mul(&temp1, S[k][i]);  
                ECP_BN462_add(&proof_theta[i][j], &temp1); 
            }
        }
   }
    for (i =0; i<2; i++){ 
        for (j =0; j<2; j++){ 
            for (k =0; k<lenX; k++){
                for (l =0; l<lenY; l++){
                    iota2(Y[l], out_iota2); 
                    ECP2_BN462_copy(&temp2,&out_iota2[j]);
                    BIG_BN462_modmul(temp, Gamma[k][l], R[k][i], r);
                    PAIR_BN462_G2mul(&temp2, temp);  
                    ECP2_BN462_add(&proof_pi[i][j], &temp2); 

                    iota1(X[k], out_iota1);       
                    ECP_BN462_copy(&temp1,&out_iota1[j]);
                    BIG_BN462_modmul(temp, Gamma[k][l], S[l][i], r);
                    PAIR_BN462_G1mul(&temp1, temp);  
                    ECP_BN462_add(&proof_theta[i][j], &temp1); 

                    BIG_BN462_modmul(temp, Gamma[k][l], S[l][j], r);
                    BIG_BN462_modmul(temp, temp, R[k][i], r) ; 
                    BIG_BN462_modadd(RGS[i][j], RGS[i][j], temp, r); 
                }
            }
        }
    }
    for (i =0; i<2; i++){ 
        for (j =0; j<2; j++){ 
            for (k =0; k<2; k++){
                ECP2_BN462_copy(&temp2,&crs_v[k][j]); 
                PAIR_BN462_G2mul(&temp2, RGS[i][k]); 
                ECP2_BN462_add(&proof_pi[i][j], &temp2);

                ECP_BN462_copy(&temp1,&crs_u[k][j]); 
                PAIR_BN462_G1mul(&temp1, T[i][k]); 
                ECP_BN462_add(&proof_theta[i][j], &temp1);
            }
        }
    }
}

bool verify_pairing_linear_XB(int lenC, ECP_BN462 crs_u[2][2], ECP2_BN462 B[], ECP_BN462 C[lenC][2], ECP2_BN462 proof[], FP12_BN462 t){ 
    int i, j, k;
    FP12_BN462 lhs, rhs, temp, out_iotat[2][2];
    ECP2_BN462 out_iota2[2];
    iotat(t,out_iotat);
    
    for( i=0; i<2 ; i++){
        for( j=0; j<2 ; j++){
            FP12_BN462_one(&lhs);
            for( k=0; k<lenC ; k++){
                iota2(B[k], out_iota2); 
                temp =PAIR_BN462_e(out_iota2[j], C[k][i]);
                FP12_BN462_ssmul(&lhs, &temp);
            }

            rhs = out_iotat[i][j];
            for( k=0; k<2 ; k++){
                iota2(proof[k], out_iota2); 
                temp =PAIR_BN462_e(out_iota2[j], crs_u[k][i]);
                FP12_BN462_ssmul(&rhs, &temp);
            }
            
            if(FP12_BN462_equals(&lhs,&rhs)==0){
                return false;
            }
        }
    }

    return true;
}

bool verify_pairing_linear_AY(int lenA, ECP2_BN462 crs_v[2][2], ECP_BN462 A[], ECP2_BN462 D[2][2], ECP_BN462 proof[], FP12_BN462 t){
    int i, j, k;
    FP12_BN462 lhs, rhs, temp, out_iotat[2][2];
    ECP_BN462 out_iota1[2];
    iotat(t,out_iotat);
    
    for( i=0; i<2 ; i++){
        for( j=0; j<2 ; j++){
            FP12_BN462_one(&lhs);
            for( k=0; k<lenA ; k++){
                iota1(A[k], out_iota1); 
                temp =PAIR_BN462_e(D[k][j],out_iota1[i]);
                FP12_BN462_ssmul(&lhs, &temp);
            }

            rhs = out_iotat[i][j];
            for( k=0; k<2 ; k++){
                iota1(proof[k], out_iota1); 
                temp =PAIR_BN462_e(crs_v[k][j],out_iota1[i]);
                FP12_BN462_ssmul(&rhs, &temp);
            }
            if(FP12_BN462_equals(&lhs,&rhs)==0){
                return false;
            }
        }
    }
    
    return true;
}

bool verify_pairing_quadratic(int lenA, int lenC, int lenD, ECP_BN462 crs_u[2][2], ECP2_BN462 crs_v[2][2], ECP_BN462 C[1][2], ECP2_BN462 D[lenD][2], ECP_BN462 A[], ECP2_BN462 B[], BIG_BN462 Gamma[1][lenD], ECP2_BN462 proof_pi[2][2], ECP_BN462 proof_theta[2][2], FP12_BN462 t){
    int i, j, k, l;
    FP12_BN462 lhs, rhs, temp, out_iotat[2][2];
    ECP2_BN462 out_iota2[2], GD[lenD], temp2; 
    ECP_BN462 out_iota1[2];
    iotat(t,out_iotat);
    
    for( i=0; i<2 ; i++){ 
        for( j=0; j<2 ; j++){ 
            FP12_BN462_one(&lhs);
            for(k=0; k<lenA ; k++){
                iota1(A[k], out_iota1); 
                temp =PAIR_BN462_e(D[k][j],out_iota1[i]);
                FP12_BN462_ssmul(&lhs, &temp);
            }

            for(k=0; k<lenC ; k++){
                iota2(B[k], out_iota2); 
                temp =PAIR_BN462_e(out_iota2[j], C[k][i]);
                FP12_BN462_ssmul(&lhs, &temp);
            }

            ECP2_BN462_inf(&GD[0]); 
            for(k=0; k<lenD ; k++){
                for(l=0; l<lenC ; l++){
                    ECP2_BN462_copy(&temp2,&D[k][j]); 
                    PAIR_BN462_G2mul(&temp2, Gamma[k][l]); 
                    ECP2_BN462_add(&GD[k], &temp2);
                }
            }

            for(k=0; k<lenC ; k++){
                temp =PAIR_BN462_e(GD[k],C[k][i]);
                FP12_BN462_ssmul(&lhs, &temp); 
            }

            rhs = out_iotat[i][j];
            for(k=0; k<2 ; k++){
                temp =PAIR_BN462_e(proof_pi[k][j],crs_u[k][i]);
                FP12_BN462_ssmul(&rhs, &temp);
                temp =PAIR_BN462_e(crs_v[k][j], proof_theta[k][i]);
                FP12_BN462_ssmul(&rhs, &temp);
            }
            if(FP12_BN462_equals(&lhs,&rhs)==0){
                return false;
            }
        }
    } 

    return true;  
}

bool verify_pairing_quadratic1(int lenA, int lenC, int lenD, ECP_BN462 crs_u[2][2], ECP2_BN462 crs_v[2][2], ECP_BN462 C[1][2], ECP2_BN462 D[2][2], ECP_BN462 A[], ECP2_BN462 B[], BIG_BN462 Gamma[1][lenD], ECP2_BN462 proof_pi[2][2], ECP_BN462 proof_theta[2][2], FP12_BN462 t){
    int i, j, k, l;
    FP12_BN462 lhs, rhs, temp, out_iotat[2][2];
    ECP2_BN462 out_iota2[2], GD[lenC], temp2; 
    ECP_BN462 out_iota1[2];
    iotat(t,out_iotat);

    for( i=0; i<2 ; i++){ 
        for( j=0; j<2 ; j++){ 
            FP12_BN462_one(&lhs);
            ECP2_BN462_inf(&GD[0]); 
            for(k=0; k<lenA ; k++){
                iota1(A[i], out_iota1); 
                temp =PAIR_BN462_e(D[k][j],out_iota1[k]);
                FP12_BN462_ssmul(&lhs, &temp);
            }

            for(k=0; k<lenC ; k++){
                iota2(B[k], out_iota2); 
                temp =PAIR_BN462_e(out_iota2[j], C[k][i]);
                FP12_BN462_ssmul(&lhs, &temp);
            }

            for(k=0; k<lenC ; k++){
                for(l=0; l<lenD ; l++){
                    ECP2_BN462_copy(&temp2,&D[k][j]); 
                    PAIR_BN462_G2mul(&temp2, Gamma[k][l]); 
                    ECP2_BN462_add(&GD[k], &temp2);
                }
                temp =PAIR_BN462_e(GD[k],C[k][i]);
                FP12_BN462_ssmul(&lhs, &temp);                
            }

            rhs = out_iotat[i][j];
            for(k=0; k<2 ; k++){
                temp =PAIR_BN462_e(proof_pi[k][j],crs_u[k][i]);
                FP12_BN462_ssmul(&rhs, &temp);
                temp =PAIR_BN462_e(crs_v[k][j], proof_theta[k][i]);
                FP12_BN462_ssmul(&rhs, &temp);
            }
               
            if(FP12_BN462_equals(&lhs,&rhs)==0){
                return false;
            }
        }
    } 
    
    return true;  
}

bool verify_pairing_quadratic2(int lenA, int lenC, int lenD, ECP_BN462 crs_u[2][2], ECP2_BN462 crs_v[2][2], ECP_BN462 C[1][2], ECP2_BN462 D[lenD][2], ECP_BN462 A[], ECP2_BN462 B[], BIG_BN462 Gamma[lenC][lenD], ECP2_BN462 proof_pi[2][2], ECP_BN462 proof_theta[2][2], FP12_BN462 t){
    int i, j, k, l;
    FP12_BN462 lhs, rhs, temp, out_iotat[2][2];
    ECP2_BN462 out_iota2[2], GD, temp2; 
    ECP_BN462 out_iota1[2];
    iotat(t,out_iotat);
    for( i=0; i<2 ; i++){ 
        for( j=0; j<2 ; j++){ 
            FP12_BN462_one(&lhs);
            ECP2_BN462_inf(&GD); 
            
            for(k=0; k<lenA ; k++){
                iota1(A[i], out_iota1); 
                temp =PAIR_BN462_e(D[k][j],out_iota1[k]);
                FP12_BN462_ssmul(&lhs, &temp);
            }

            for(k=0; k<lenC ; k++){
                iota2(B[k], out_iota2); 
                temp =PAIR_BN462_e(out_iota2[j], C[k][i]);
                FP12_BN462_ssmul(&lhs, &temp);
            }

            for(k=0; k<lenC ; k++){
                for(l=0; l<lenD ; l++){
                    ECP2_BN462_copy(&temp2,&D[l][j]); 
                    PAIR_BN462_G2mul(&temp2, Gamma[k][l]); 
                    ECP2_BN462_copy(&GD, &temp2);
                }
                temp =PAIR_BN462_e(GD,C[k][i]); 
                FP12_BN462_ssmul(&lhs, &temp);                
            }

            rhs = out_iotat[i][j];
            for(k=0; k<2 ; k++){
                temp =PAIR_BN462_e(proof_pi[k][j],crs_u[k][i]);
                FP12_BN462_ssmul(&rhs, &temp);
                temp =PAIR_BN462_e(crs_v[k][j], proof_theta[k][i]);
                FP12_BN462_ssmul(&rhs, &temp);
            }
            
            if(FP12_BN462_equals(&lhs,&rhs)==0){
                return false;
            }
        }
    }

    return true;  
}
