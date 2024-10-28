#include "fields_cos.h"
#include <math.h>

double Sigma(double *ypar, double *ycos, double *ysin, double L) {
   double Sigma_result;
   Sigma_result = ycos[3]*ypar[15] + ycos[4]*ypar[15] + ycos[5]*ypar[15] + ycos[3]*ypar[16]*ycos[10] + ycos[3]*ypar[16]*ycos[11] + ycos[4]*ypar[16]*ycos[11] + ycos[5]*ypar[16]*ycos[10] + ycos[9]*ycos[4]*ypar[16] + ycos[9]*ycos[5]*ypar[16] + ycos[3]*ycos[5]*ycos[4]*ypar[17] + ycos[3]*ypar[18]*ycos[10]*ycos[11] + ycos[9]*ycos[4]*ypar[18]*ycos[11] + ycos[9]*ycos[5]*ypar[18]*ycos[10];
   return Sigma_result;
}

double Pi1(double *ypar, double *ycos, double *ysin, double L) {
   double Pi1_result;
   Pi1_result = -ycos[7]*ypar[10]*ysin[12] + ypar[13]*ysin[12]*ycos[13] + ypar[5]*ycos[1]*ysin[12] + ysin[0]*ypar[0]*ycos[1] + ysin[0]*ypar[2]*ycos[7] + ysin[0]*ypar[5]*ycos[13] + ysin[6]*ypar[10]*ycos[13] - ysin[6]*ypar[2]*ycos[1] + ysin[6]*ypar[8]*ycos[7] + ycos[7]*ycos[5]*ypar[11]*ysin[12] - ycos[7]*ypar[12]*ysin[12]*ycos[11] + ypar[14]*ysin[12]*ycos[11]*ycos[13] - ypar[6]*ycos[1]*ycos[5]*ysin[12] + ypar[7]*ycos[1]*ysin[12]*ycos[11] + ysin[0]*ypar[1]*ycos[1]*ycos[11] + ysin[0]*ypar[3]*ycos[7]*ycos[5] + ysin[0]*ypar[4]*ycos[7]*ycos[11] + ysin[0]*ypar[6]*ycos[5]*ycos[13] + ysin[0]*ypar[7]*ycos[11]*ycos[13] + ysin[6]*ycos[5]*ypar[11]*ycos[13] + ysin[6]*ypar[12]*ycos[11]*ycos[13] + ysin[6]*ypar[3]*ycos[1]*ycos[5] - ysin[6]*ypar[4]*ycos[1]*ycos[11] + ysin[6]*ypar[9]*ycos[7]*ycos[11];
   return Pi1_result;
}

double Pi2(double *ypar, double *ycos, double *ysin, double L) {
   double Pi2_result;
   Pi2_result = -ycos[8]*ypar[10]*ysin[13] + ypar[13]*ysin[13]*ycos[14] + ypar[5]*ycos[2]*ysin[13] + ysin[1]*ypar[0]*ycos[2] + ysin[1]*ypar[2]*ycos[8] + ysin[1]*ypar[5]*ycos[14] + ysin[7]*ypar[10]*ycos[14] - ysin[7]*ypar[2]*ycos[2] + ysin[7]*ypar[8]*ycos[8] + ycos[3]*ycos[8]*ypar[11]*ysin[13] - ycos[9]*ycos[8]*ypar[12]*ysin[13] + ycos[9]*ypar[14]*ysin[13]*ycos[14] - ypar[6]*ycos[3]*ycos[2]*ysin[13] + ypar[7]*ycos[2]*ycos[9]*ysin[13] + ysin[1]*ypar[1]*ycos[2]*ycos[9] + ysin[1]*ypar[3]*ycos[3]*ycos[8] + ysin[1]*ypar[4]*ycos[9]*ycos[8] + ysin[1]*ypar[6]*ycos[3]*ycos[14] + ysin[1]*ypar[7]*ycos[9]*ycos[14] + ysin[7]*ycos[3]*ypar[11]*ycos[14] + ysin[7]*ycos[9]*ypar[12]*ycos[14] + ysin[7]*ypar[3]*ycos[3]*ycos[2] - ysin[7]*ypar[4]*ycos[2]*ycos[9] + ysin[7]*ypar[9]*ycos[9]*ycos[8];
   return Pi2_result;
}

double Pi3(double *ypar, double *ycos, double *ysin, double L) {
   double Pi3_result;
   Pi3_result = -ycos[6]*ypar[10]*ysin[14] + ypar[13]*ysin[14]*ycos[12] + ypar[5]*ycos[0]*ysin[14] + ysin[2]*ypar[0]*ycos[0] + ysin[2]*ypar[2]*ycos[6] + ysin[2]*ypar[5]*ycos[12] + ysin[8]*ypar[10]*ycos[12] - ysin[8]*ypar[2]*ycos[0] + ysin[8]*ypar[8]*ycos[6] + ycos[6]*ycos[4]*ypar[11]*ysin[14] - ycos[6]*ypar[12]*ysin[14]*ycos[10] + ypar[14]*ysin[14]*ycos[10]*ycos[12] - ypar[6]*ycos[0]*ycos[4]*ysin[14] + ypar[7]*ycos[0]*ysin[14]*ycos[10] + ysin[2]*ypar[1]*ycos[0]*ycos[10] + ysin[2]*ypar[3]*ycos[6]*ycos[4] + ysin[2]*ypar[4]*ycos[6]*ycos[10] + ysin[2]*ypar[6]*ycos[4]*ycos[12] + ysin[2]*ypar[7]*ycos[10]*ycos[12] + ysin[8]*ycos[4]*ypar[11]*ycos[12] + ysin[8]*ypar[12]*ycos[10]*ycos[12] + ysin[8]*ypar[3]*ycos[0]*ycos[4] - ysin[8]*ypar[4]*ycos[0]*ycos[10] + ysin[8]*ypar[9]*ycos[6]*ycos[10];
   return Pi3_result;
}

double Sigma_1(double *ypar, double *ycos, double *ysin, double L) {
   double Sigma_1_result;
   Sigma_1_result = -2*M_PI*ysin[3]*ypar[15]/L - 2*M_PI*ysin[3]*ypar[16]*ycos[10]/L - 2*M_PI*ysin[3]*ypar[16]*ycos[11]/L - 4*M_PI*ysin[9]*ycos[4]*ypar[16]/L - 4*M_PI*ysin[9]*ycos[5]*ypar[16]/L - 2*M_PI*ysin[3]*ycos[5]*ycos[4]*ypar[17]/L - 2*M_PI*ysin[3]*ypar[18]*ycos[10]*ycos[11]/L - 4*M_PI*ysin[9]*ycos[4]*ypar[18]*ycos[11]/L - 4*M_PI*ysin[9]*ycos[5]*ypar[18]*ycos[10]/L;
   return Sigma_1_result;
}

double Sigma_2(double *ypar, double *ycos, double *ysin, double L) {
   double Sigma_2_result;
   Sigma_2_result = -2*M_PI*ysin[4]*ypar[15]/L - 4*M_PI*ycos[3]*ypar[16]*ysin[10]/L - 4*M_PI*ycos[5]*ypar[16]*ysin[10]/L - 2*M_PI*ysin[4]*ycos[9]*ypar[16]/L - 2*M_PI*ysin[4]*ypar[16]*ycos[11]/L - 4*M_PI*ycos[3]*ypar[18]*ysin[10]*ycos[11]/L - 4*M_PI*ycos[9]*ycos[5]*ypar[18]*ysin[10]/L - 2*M_PI*ysin[4]*ycos[3]*ycos[5]*ypar[17]/L - 2*M_PI*ysin[4]*ycos[9]*ypar[18]*ycos[11]/L;
   return Sigma_2_result;
}

double Sigma_3(double *ypar, double *ycos, double *ysin, double L) {
   double Sigma_3_result;
   Sigma_3_result = -2*M_PI*ysin[5]*ypar[15]/L - 4*M_PI*ycos[3]*ypar[16]*ysin[11]/L - 4*M_PI*ycos[4]*ypar[16]*ysin[11]/L - 2*M_PI*ysin[5]*ycos[9]*ypar[16]/L - 2*M_PI*ysin[5]*ypar[16]*ycos[10]/L - 4*M_PI*ycos[3]*ypar[18]*ysin[11]*ycos[10]/L - 4*M_PI*ycos[9]*ycos[4]*ypar[18]*ysin[11]/L - 2*M_PI*ysin[5]*ycos[3]*ycos[4]*ypar[17]/L - 2*M_PI*ysin[5]*ycos[9]*ypar[18]*ycos[10]/L;
   return Sigma_3_result;
}

double Pi1_1(double *ypar, double *ycos, double *ysin, double L) {
   double Pi1_1_result;
   Pi1_1_result = 3*M_PI*ycos[6]*ypar[10]*ycos[13]/L - 5*M_PI*ycos[7]*ypar[10]*ycos[12]/L + M_PI*ypar[0]*ycos[1]*ycos[0]/L + 5*M_PI*ypar[13]*ycos[13]*ycos[12]/L + M_PI*ypar[2]*ycos[0]*ycos[7]/L - 3*M_PI*ypar[2]*ycos[1]*ycos[6]/L + M_PI*ypar[5]*ycos[0]*ycos[13]/L + 5*M_PI*ypar[5]*ycos[1]*ycos[12]/L + 3*M_PI*ypar[8]*ycos[7]*ycos[6]/L + 3*M_PI*ycos[6]*ycos[5]*ypar[11]*ycos[13]/L + 3*M_PI*ycos[6]*ypar[12]*ycos[11]*ycos[13]/L + 5*M_PI*ycos[7]*ycos[5]*ypar[11]*ycos[12]/L - 5*M_PI*ycos[7]*ypar[12]*ycos[11]*ycos[12]/L + 5*M_PI*ypar[14]*ycos[11]*ycos[13]*ycos[12]/L + M_PI*ypar[1]*ycos[1]*ycos[0]*ycos[11]/L + M_PI*ypar[3]*ycos[0]*ycos[7]*ycos[5]/L + 3*M_PI*ypar[3]*ycos[1]*ycos[6]*ycos[5]/L + M_PI*ypar[4]*ycos[0]*ycos[7]*ycos[11]/L - 3*M_PI*ypar[4]*ycos[1]*ycos[6]*ycos[11]/L + M_PI*ypar[6]*ycos[0]*ycos[5]*ycos[13]/L - 5*M_PI*ypar[6]*ycos[1]*ycos[5]*ycos[12]/L + M_PI*ypar[7]*ycos[0]*ycos[11]*ycos[13]/L + 5*M_PI*ypar[7]*ycos[1]*ycos[11]*ycos[12]/L + 3*M_PI*ypar[9]*ycos[7]*ycos[6]*ycos[11]/L;
   return Pi1_1_result;
}

double Pi1_2(double *ypar,double *ycos, double *ysin, double L) {
   double Pi1_2_result;
   Pi1_2_result = -5*M_PI*ypar[13]*ysin[13]*ysin[12]/L - 5*M_PI*ysin[0]*ypar[5]*ysin[13]/L - M_PI*ysin[0]*ysin[1]*ypar[0]/L - 3*M_PI*ysin[0]*ysin[7]*ypar[2]/L - M_PI*ysin[1]*ypar[5]*ysin[12]/L + M_PI*ysin[1]*ysin[6]*ypar[2]/L - 5*M_PI*ysin[6]*ypar[10]*ysin[13]/L - 3*M_PI*ysin[6]*ysin[7]*ypar[8]/L + 3*M_PI*ysin[7]*ypar[10]*ysin[12]/L - 5*M_PI*ypar[14]*ysin[13]*ysin[12]*ycos[11]/L - 5*M_PI*ysin[0]*ypar[6]*ycos[5]*ysin[13]/L - 5*M_PI*ysin[0]*ypar[7]*ysin[13]*ycos[11]/L - M_PI*ysin[0]*ysin[1]*ypar[1]*ycos[11]/L - 3*M_PI*ysin[0]*ysin[7]*ypar[3]*ycos[5]/L - 3*M_PI*ysin[0]*ysin[7]*ypar[4]*ycos[11]/L + M_PI*ysin[1]*ypar[6]*ycos[5]*ysin[12]/L - M_PI*ysin[1]*ypar[7]*ysin[12]*ycos[11]/L - M_PI*ysin[1]*ysin[6]*ypar[3]*ycos[5]/L + M_PI*ysin[1]*ysin[6]*ypar[4]*ycos[11]/L - 5*M_PI*ysin[6]*ycos[5]*ypar[11]*ysin[13]/L - 5*M_PI*ysin[6]*ypar[12]*ysin[13]*ycos[11]/L - 3*M_PI*ysin[6]*ysin[7]*ypar[9]*ycos[11]/L - 3*M_PI*ysin[7]*ycos[5]*ypar[11]*ysin[12]/L + 3*M_PI*ysin[7]*ypar[12]*ysin[12]*ycos[11]/L;
   return Pi1_2_result;
}

double Pi1_3(double *ypar,double *ycos, double *ysin, double L) {
   double Pi1_3_result;
   Pi1_3_result = 4*M_PI*ycos[7]*ypar[12]*ysin[11]*ysin[12]/L - 4*M_PI*ypar[14]*ysin[11]*ysin[12]*ycos[13]/L - 4*M_PI*ypar[7]*ycos[1]*ysin[11]*ysin[12]/L - 4*M_PI*ysin[0]*ypar[1]*ycos[1]*ysin[11]/L - 4*M_PI*ysin[0]*ypar[4]*ycos[7]*ysin[11]/L - 4*M_PI*ysin[0]*ypar[7]*ysin[11]*ycos[13]/L - 2*M_PI*ysin[0]*ysin[5]*ypar[3]*ycos[7]/L - 2*M_PI*ysin[0]*ysin[5]*ypar[6]*ycos[13]/L - 2*M_PI*ysin[5]*ycos[7]*ypar[11]*ysin[12]/L + 2*M_PI*ysin[5]*ypar[6]*ycos[1]*ysin[12]/L - 4*M_PI*ysin[6]*ypar[12]*ysin[11]*ycos[13]/L + 4*M_PI*ysin[6]*ypar[4]*ycos[1]*ysin[11]/L - 4*M_PI*ysin[6]*ypar[9]*ycos[7]*ysin[11]/L - 2*M_PI*ysin[6]*ysin[5]*ypar[11]*ycos[13]/L - 2*M_PI*ysin[6]*ysin[5]*ypar[3]*ycos[1]/L;
   return Pi1_3_result;
}

double Pi2_1(double *ypar,double *ycos, double *ysin, double L) {
   double Pi2_1_result;
   Pi2_1_result = -2*M_PI*ysin[1]*ysin[3]*ypar[3]*ycos[8]/L - 2*M_PI*ysin[1]*ysin[3]*ypar[6]*ycos[14]/L - 4*M_PI*ysin[1]*ysin[9]*ypar[1]*ycos[2]/L - 4*M_PI*ysin[1]*ysin[9]*ypar[4]*ycos[8]/L - 4*M_PI*ysin[1]*ysin[9]*ypar[7]*ycos[14]/L - 2*M_PI*ysin[3]*ycos[8]*ypar[11]*ysin[13]/L + 2*M_PI*ysin[3]*ypar[6]*ycos[2]*ysin[13]/L - 2*M_PI*ysin[3]*ysin[7]*ypar[11]*ycos[14]/L - 2*M_PI*ysin[3]*ysin[7]*ypar[3]*ycos[2]/L - 4*M_PI*ysin[7]*ysin[9]*ypar[12]*ycos[14]/L + 4*M_PI*ysin[7]*ysin[9]*ypar[4]*ycos[2]/L - 4*M_PI*ysin[7]*ysin[9]*ypar[9]*ycos[8]/L + 4*M_PI*ysin[9]*ycos[8]*ypar[12]*ysin[13]/L - 4*M_PI*ysin[9]*ypar[14]*ysin[13]*ycos[14]/L - 4*M_PI*ysin[9]*ypar[7]*ycos[2]*ysin[13]/L;
   return Pi2_1_result;
}

double Pi2_2(double *ypar, double *ycos, double *ysin, double L) {
   double Pi2_2_result;
   Pi2_2_result = 3*M_PI*ycos[7]*ypar[10]*ycos[14]/L - 5*M_PI*ycos[8]*ypar[10]*ycos[13]/L + M_PI*ypar[0]*ycos[1]*ycos[2]/L + 5*M_PI*ypar[13]*ycos[13]*ycos[14]/L + M_PI*ypar[2]*ycos[1]*ycos[8]/L - 3*M_PI*ypar[2]*ycos[2]*ycos[7]/L + M_PI*ypar[5]*ycos[1]*ycos[14]/L + 5*M_PI*ypar[5]*ycos[2]*ycos[13]/L + 3*M_PI*ypar[8]*ycos[7]*ycos[8]/L + 3*M_PI*ycos[3]*ycos[7]*ypar[11]*ycos[14]/L + 5*M_PI*ycos[3]*ycos[8]*ypar[11]*ycos[13]/L + 3*M_PI*ycos[7]*ycos[9]*ypar[12]*ycos[14]/L - 5*M_PI*ycos[9]*ycos[8]*ypar[12]*ycos[13]/L + 5*M_PI*ycos[9]*ypar[14]*ycos[13]*ycos[14]/L + M_PI*ypar[1]*ycos[1]*ycos[2]*ycos[9]/L + M_PI*ypar[3]*ycos[1]*ycos[3]*ycos[8]/L + 3*M_PI*ypar[3]*ycos[3]*ycos[2]*ycos[7]/L + M_PI*ypar[4]*ycos[1]*ycos[9]*ycos[8]/L - 3*M_PI*ypar[4]*ycos[2]*ycos[7]*ycos[9]/L + M_PI*ypar[6]*ycos[1]*ycos[3]*ycos[14]/L - 5*M_PI*ypar[6]*ycos[3]*ycos[2]*ycos[13]/L + M_PI*ypar[7]*ycos[1]*ycos[9]*ycos[14]/L + 5*M_PI*ypar[7]*ycos[2]*ycos[9]*ycos[13]/L + 3*M_PI*ypar[9]*ycos[7]*ycos[9]*ycos[8]/L;
   return Pi2_2_result;
}

double Pi2_3(double *ypar,double *ycos, double *ysin, double L) {
   double Pi2_3_result;
   Pi2_3_result = -5*M_PI*ypar[13]*ysin[13]*ysin[14]/L - 5*M_PI*ysin[1]*ypar[5]*ysin[14]/L - 3*M_PI*ysin[1]*ysin[8]*ypar[2]/L - M_PI*ysin[2]*ypar[5]*ysin[13]/L - M_PI*ysin[2]*ysin[1]*ypar[0]/L + M_PI*ysin[2]*ysin[7]*ypar[2]/L - 5*M_PI*ysin[7]*ypar[10]*ysin[14]/L + 3*M_PI*ysin[8]*ypar[10]*ysin[13]/L - 3*M_PI*ysin[8]*ysin[7]*ypar[8]/L - 5*M_PI*ycos[9]*ypar[14]*ysin[13]*ysin[14]/L - 5*M_PI*ysin[1]*ypar[6]*ycos[3]*ysin[14]/L - 5*M_PI*ysin[1]*ypar[7]*ycos[9]*ysin[14]/L - 3*M_PI*ysin[1]*ysin[8]*ypar[3]*ycos[3]/L - 3*M_PI*ysin[1]*ysin[8]*ypar[4]*ycos[9]/L + M_PI*ysin[2]*ypar[6]*ycos[3]*ysin[13]/L - M_PI*ysin[2]*ypar[7]*ycos[9]*ysin[13]/L - M_PI*ysin[2]*ysin[1]*ypar[1]*ycos[9]/L - M_PI*ysin[2]*ysin[7]*ypar[3]*ycos[3]/L + M_PI*ysin[2]*ysin[7]*ypar[4]*ycos[9]/L - 5*M_PI*ysin[7]*ycos[3]*ypar[11]*ysin[14]/L - 5*M_PI*ysin[7]*ycos[9]*ypar[12]*ysin[14]/L - 3*M_PI*ysin[8]*ycos[3]*ypar[11]*ysin[13]/L + 3*M_PI*ysin[8]*ycos[9]*ypar[12]*ysin[13]/L - 3*M_PI*ysin[8]*ysin[7]*ypar[9]*ycos[9]/L;
   return Pi2_3_result;
}

double Pi3_1(double *ypar,double *ycos, double *ysin, double L) {
   double Pi3_1_result;
   Pi3_1_result = -5*M_PI*ypar[13]*ysin[14]*ysin[12]/L - M_PI*ysin[0]*ypar[5]*ysin[14]/L - M_PI*ysin[0]*ysin[2]*ypar[0]/L + M_PI*ysin[0]*ysin[8]*ypar[2]/L - 5*M_PI*ysin[2]*ypar[5]*ysin[12]/L - 3*M_PI*ysin[2]*ysin[6]*ypar[2]/L + 3*M_PI*ysin[6]*ypar[10]*ysin[14]/L - 3*M_PI*ysin[6]*ysin[8]*ypar[8]/L - 5*M_PI*ysin[8]*ypar[10]*ysin[12]/L - 5*M_PI*ypar[14]*ysin[14]*ysin[12]*ycos[10]/L + M_PI*ysin[0]*ypar[6]*ycos[4]*ysin[14]/L - M_PI*ysin[0]*ypar[7]*ysin[14]*ycos[10]/L - M_PI*ysin[0]*ysin[2]*ypar[1]*ycos[10]/L - M_PI*ysin[0]*ysin[8]*ypar[3]*ycos[4]/L + M_PI*ysin[0]*ysin[8]*ypar[4]*ycos[10]/L - 5*M_PI*ysin[2]*ypar[6]*ycos[4]*ysin[12]/L - 5*M_PI*ysin[2]*ypar[7]*ysin[12]*ycos[10]/L - 3*M_PI*ysin[2]*ysin[6]*ypar[3]*ycos[4]/L - 3*M_PI*ysin[2]*ysin[6]*ypar[4]*ycos[10]/L - 3*M_PI*ysin[6]*ycos[4]*ypar[11]*ysin[14]/L + 3*M_PI*ysin[6]*ypar[12]*ysin[14]*ycos[10]/L - 3*M_PI*ysin[6]*ysin[8]*ypar[9]*ycos[10]/L - 5*M_PI*ysin[8]*ycos[4]*ypar[11]*ysin[12]/L - 5*M_PI*ysin[8]*ypar[12]*ysin[12]*ycos[10]/L;
   return Pi3_1_result;
}

double Pi3_2(double *ypar,double *ycos, double *ysin, double L) {
   double Pi3_2_result;
   Pi3_2_result = 4*M_PI*ycos[6]*ypar[12]*ysin[14]*ysin[10]/L - 4*M_PI*ypar[14]*ysin[14]*ysin[10]*ycos[12]/L - 4*M_PI*ypar[7]*ycos[0]*ysin[14]*ysin[10]/L - 4*M_PI*ysin[2]*ypar[1]*ycos[0]*ysin[10]/L - 4*M_PI*ysin[2]*ypar[4]*ycos[6]*ysin[10]/L - 4*M_PI*ysin[2]*ypar[7]*ysin[10]*ycos[12]/L - 2*M_PI*ysin[2]*ysin[4]*ypar[3]*ycos[6]/L - 2*M_PI*ysin[2]*ysin[4]*ypar[6]*ycos[12]/L - 2*M_PI*ysin[4]*ycos[6]*ypar[11]*ysin[14]/L + 2*M_PI*ysin[4]*ypar[6]*ycos[0]*ysin[14]/L - 2*M_PI*ysin[4]*ysin[8]*ypar[11]*ycos[12]/L - 2*M_PI*ysin[4]*ysin[8]*ypar[3]*ycos[0]/L - 4*M_PI*ysin[8]*ypar[12]*ysin[10]*ycos[12]/L + 4*M_PI*ysin[8]*ypar[4]*ycos[0]*ysin[10]/L - 4*M_PI*ysin[8]*ypar[9]*ycos[6]*ysin[10]/L;
   return Pi3_2_result;
}

double Pi3_3(double *ypar, double *ycos, double *ysin, double L) {
   double Pi3_3_result;
   Pi3_3_result = -5*M_PI*ycos[6]*ypar[10]*ycos[14]/L + 3*M_PI*ycos[8]*ypar[10]*ycos[12]/L + M_PI*ypar[0]*ycos[0]*ycos[2]/L + 5*M_PI*ypar[13]*ycos[14]*ycos[12]/L - 3*M_PI*ypar[2]*ycos[0]*ycos[8]/L + M_PI*ypar[2]*ycos[2]*ycos[6]/L + 5*M_PI*ypar[5]*ycos[0]*ycos[14]/L + M_PI*ypar[5]*ycos[2]*ycos[12]/L + 3*M_PI*ypar[8]*ycos[6]*ycos[8]/L + 5*M_PI*ycos[6]*ycos[4]*ypar[11]*ycos[14]/L - 5*M_PI*ycos[6]*ypar[12]*ycos[10]*ycos[14]/L + 3*M_PI*ycos[8]*ycos[4]*ypar[11]*ycos[12]/L + 3*M_PI*ycos[8]*ypar[12]*ycos[10]*ycos[12]/L + 5*M_PI*ypar[14]*ycos[10]*ycos[14]*ycos[12]/L + M_PI*ypar[1]*ycos[0]*ycos[2]*ycos[10]/L + 3*M_PI*ypar[3]*ycos[0]*ycos[8]*ycos[4]/L + M_PI*ypar[3]*ycos[2]*ycos[6]*ycos[4]/L - 3*M_PI*ypar[4]*ycos[0]*ycos[8]*ycos[10]/L + M_PI*ypar[4]*ycos[2]*ycos[6]*ycos[10]/L - 5*M_PI*ypar[6]*ycos[0]*ycos[4]*ycos[14]/L + M_PI*ypar[6]*ycos[2]*ycos[4]*ycos[12]/L + 5*M_PI*ypar[7]*ycos[0]*ycos[10]*ycos[14]/L + M_PI*ypar[7]*ycos[2]*ycos[10]*ycos[12]/L + 3*M_PI*ypar[9]*ycos[6]*ycos[8]*ycos[10]/L;
   return Pi3_3_result;
}

