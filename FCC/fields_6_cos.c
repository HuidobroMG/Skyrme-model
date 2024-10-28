#include "fields_cos.h"
#include <math.h>

double Sigma(double *ypar, double *ycos, double *ysin, double L) {
   double Sigma_result;
   Sigma_result = ycos[3]*ypar[19] + ycos[4]*ypar[13] + ycos[5]*ypar[12] + ycos[3]*ycos[4]*ypar[21] + ycos[3]*ycos[5]*ypar[20] + ycos[5]*ycos[4]*ypar[14] + ycos[0]*ycos[2]*ycos[7]*ypar[17] + ycos[0]*ycos[7]*ycos[8]*ypar[18] + ycos[1]*ycos[0]*ycos[2]*ypar[15] + ycos[1]*ycos[0]*ycos[8]*ypar[16] + ycos[1]*ycos[2]*ycos[6]*ypar[23] + ycos[1]*ycos[6]*ycos[8]*ypar[24] + ycos[2]*ycos[7]*ycos[6]*ypar[25] + ycos[3]*ycos[5]*ycos[4]*ypar[22] + ycos[7]*ycos[6]*ycos[8]*ypar[26] - (ypar[12] + ypar[13] + ypar[14] + ypar[19] + ypar[20] + ypar[21] + ypar[22]);
   return Sigma_result;
}

double Pi1(double *ypar, double *ycos, double *ysin, double L) {
   double Pi1_result;
   Pi1_result = ysin[0]*ypar[0] + ysin[6]*ypar[8] + ysin[0]*ypar[1]*ycos[5] + ysin[0]*ypar[2]*ycos[4] + ysin[6]*ycos[4]*ypar[10] + ysin[6]*ypar[9]*ycos[5] + ysin[0]*ypar[3]*ycos[5]*ycos[4] + ysin[3]*ypar[4]*ycos[1]*ycos[2] + ysin[3]*ypar[5]*ycos[1]*ycos[8] + ysin[3]*ypar[6]*ycos[2]*ycos[7] + ysin[3]*ypar[7]*ycos[7]*ycos[8] + ysin[6]*ycos[5]*ycos[4]*ypar[11];
   return Pi1_result;
}

double Pi2(double *ypar, double *ycos, double *ysin, double L) {
   double Pi2_result;
   Pi2_result = ysin[1]*ypar[0] + ysin[7]*ypar[8] + ysin[1]*ypar[1]*ycos[3] + ysin[1]*ypar[2]*ycos[5] + ysin[7]*ycos[5]*ypar[10] + ysin[7]*ypar[9]*ycos[3] + ysin[1]*ypar[3]*ycos[3]*ycos[5] + ysin[4]*ypar[4]*ycos[0]*ycos[2] + ysin[4]*ypar[5]*ycos[2]*ycos[6] + ysin[4]*ypar[6]*ycos[0]*ycos[8] + ysin[4]*ypar[7]*ycos[6]*ycos[8] + ysin[7]*ycos[3]*ycos[5]*ypar[11];
   return Pi2_result;
}

double Pi3(double *ypar, double *ycos, double *ysin, double L) {
   double Pi3_result;
   Pi3_result = ysin[2]*ypar[0] + ysin[8]*ypar[8] + ysin[2]*ypar[1]*ycos[4] + ysin[2]*ypar[2]*ycos[3] + ysin[8]*ycos[3]*ypar[10] + ysin[8]*ypar[9]*ycos[4] + ysin[2]*ypar[3]*ycos[3]*ycos[4] + ysin[5]*ypar[4]*ycos[1]*ycos[0] + ysin[5]*ypar[5]*ycos[0]*ycos[7] + ysin[5]*ypar[6]*ycos[1]*ycos[6] + ysin[5]*ypar[7]*ycos[7]*ycos[6] + ysin[8]*ycos[3]*ycos[4]*ypar[11];
   return Pi3_result;
}

double Sigma_1(double *ypar, double *ycos, double *ysin, double L) {
   double Sigma_1_result;
   Sigma_1_result = -2*M_PI*ysin[3]*ypar[19]/L - 2*M_PI*ysin[3]*ycos[4]*ypar[21]/L - 2*M_PI*ysin[3]*ycos[5]*ypar[20]/L - M_PI*ysin[0]*ycos[1]*ycos[2]*ypar[15]/L - M_PI*ysin[0]*ycos[1]*ycos[8]*ypar[16]/L - M_PI*ysin[0]*ycos[2]*ycos[7]*ypar[17]/L - M_PI*ysin[0]*ycos[7]*ycos[8]*ypar[18]/L - 2*M_PI*ysin[3]*ycos[5]*ycos[4]*ypar[22]/L - 3*M_PI*ysin[6]*ycos[1]*ycos[2]*ypar[23]/L - 3*M_PI*ysin[6]*ycos[1]*ycos[8]*ypar[24]/L - 3*M_PI*ysin[6]*ycos[2]*ycos[7]*ypar[25]/L - 3*M_PI*ysin[6]*ycos[7]*ycos[8]*ypar[26]/L;
   return Sigma_1_result;
}

double Sigma_2(double *ypar, double *ycos, double *ysin, double L) {
   double Sigma_2_result;
   Sigma_2_result = -2*M_PI*ysin[4]*ypar[13]/L - 2*M_PI*ysin[4]*ycos[3]*ypar[21]/L - 2*M_PI*ysin[4]*ycos[5]*ypar[14]/L - M_PI*ysin[1]*ycos[0]*ycos[2]*ypar[15]/L - M_PI*ysin[1]*ycos[0]*ycos[8]*ypar[16]/L - M_PI*ysin[1]*ycos[2]*ycos[6]*ypar[23]/L - M_PI*ysin[1]*ycos[6]*ycos[8]*ypar[24]/L - 2*M_PI*ysin[4]*ycos[3]*ycos[5]*ypar[22]/L - 3*M_PI*ysin[7]*ycos[0]*ycos[2]*ypar[17]/L - 3*M_PI*ysin[7]*ycos[0]*ycos[8]*ypar[18]/L - 3*M_PI*ysin[7]*ycos[2]*ycos[6]*ypar[25]/L - 3*M_PI*ysin[7]*ycos[6]*ycos[8]*ypar[26]/L;
   return Sigma_2_result;
}

double Sigma_3(double *ypar, double *ycos, double *ysin, double L) {
   double Sigma_3_result;
   Sigma_3_result = -2*M_PI*ysin[5]*ypar[12]/L - 2*M_PI*ysin[5]*ycos[3]*ypar[20]/L - 2*M_PI*ysin[5]*ycos[4]*ypar[14]/L - M_PI*ysin[2]*ycos[0]*ycos[7]*ypar[17]/L - M_PI*ysin[2]*ycos[1]*ycos[0]*ypar[15]/L - M_PI*ysin[2]*ycos[1]*ycos[6]*ypar[23]/L - M_PI*ysin[2]*ycos[7]*ycos[6]*ypar[25]/L - 2*M_PI*ysin[5]*ycos[3]*ycos[4]*ypar[22]/L - 3*M_PI*ysin[8]*ycos[0]*ycos[7]*ypar[18]/L - 3*M_PI*ysin[8]*ycos[1]*ycos[0]*ypar[16]/L - 3*M_PI*ysin[8]*ycos[1]*ycos[6]*ypar[24]/L - 3*M_PI*ysin[8]*ycos[7]*ycos[6]*ypar[26]/L;
   return Sigma_3_result;
}

double Pi1_1(double *ypar, double *ycos, double *ysin, double L) {
   double Pi1_1_result;
   Pi1_1_result = M_PI*ypar[0]*ycos[0]/L + 3*M_PI*ypar[8]*ycos[6]/L + 3*M_PI*ycos[6]*ycos[4]*ypar[10]/L + M_PI*ypar[1]*ycos[0]*ycos[5]/L + M_PI*ypar[2]*ycos[0]*ycos[4]/L + 3*M_PI*ypar[9]*ycos[6]*ycos[5]/L + 3*M_PI*ycos[6]*ycos[5]*ycos[4]*ypar[11]/L + M_PI*ypar[3]*ycos[0]*ycos[5]*ycos[4]/L + 2*M_PI*ypar[4]*ycos[1]*ycos[3]*ycos[2]/L + 2*M_PI*ypar[5]*ycos[1]*ycos[3]*ycos[8]/L + 2*M_PI*ypar[6]*ycos[3]*ycos[2]*ycos[7]/L + 2*M_PI*ypar[7]*ycos[3]*ycos[7]*ycos[8]/L;
   return Pi1_1_result;
}

double Pi1_2(double *ypar,double *ycos, double *ysin, double L) {
   double Pi1_2_result;
   Pi1_2_result = -2*M_PI*ysin[0]*ysin[4]*ypar[2]/L - 2*M_PI*ysin[4]*ysin[6]*ypar[10]/L - 2*M_PI*ysin[0]*ysin[4]*ypar[3]*ycos[5]/L - M_PI*ysin[1]*ysin[3]*ypar[4]*ycos[2]/L - M_PI*ysin[1]*ysin[3]*ypar[5]*ycos[8]/L - 3*M_PI*ysin[3]*ysin[7]*ypar[6]*ycos[2]/L - 3*M_PI*ysin[3]*ysin[7]*ypar[7]*ycos[8]/L - 2*M_PI*ysin[4]*ysin[6]*ycos[5]*ypar[11]/L;
   return Pi1_2_result;
}

double Pi1_3(double *ypar,double *ycos, double *ysin, double L) {
   double Pi1_3_result;
   Pi1_3_result = -2*M_PI*ysin[0]*ysin[5]*ypar[1]/L - 2*M_PI*ysin[6]*ysin[5]*ypar[9]/L - 2*M_PI*ysin[0]*ysin[5]*ypar[3]*ycos[4]/L - M_PI*ysin[2]*ysin[3]*ypar[4]*ycos[1]/L - M_PI*ysin[2]*ysin[3]*ypar[6]*ycos[7]/L - 3*M_PI*ysin[3]*ysin[8]*ypar[5]*ycos[1]/L - 3*M_PI*ysin[3]*ysin[8]*ypar[7]*ycos[7]/L - 2*M_PI*ysin[6]*ysin[5]*ycos[4]*ypar[11]/L;
   return Pi1_3_result;
}

double Pi2_1(double *ypar,double *ycos, double *ysin, double L) {
   double Pi2_1_result;
   Pi2_1_result = -2*M_PI*ysin[1]*ysin[3]*ypar[1]/L - 2*M_PI*ysin[3]*ysin[7]*ypar[9]/L - M_PI*ysin[0]*ysin[4]*ypar[4]*ycos[2]/L - M_PI*ysin[0]*ysin[4]*ypar[6]*ycos[8]/L - 2*M_PI*ysin[1]*ysin[3]*ypar[3]*ycos[5]/L - 2*M_PI*ysin[3]*ysin[7]*ycos[5]*ypar[11]/L - 3*M_PI*ysin[4]*ysin[6]*ypar[5]*ycos[2]/L - 3*M_PI*ysin[4]*ysin[6]*ypar[7]*ycos[8]/L;
   return Pi2_1_result;
}

double Pi2_2(double *ypar, double *ycos, double *ysin, double L) {
   double Pi2_2_result;
   Pi2_2_result = M_PI*ypar[0]*ycos[1]/L + 3*M_PI*ypar[8]*ycos[7]/L + 3*M_PI*ycos[7]*ycos[5]*ypar[10]/L + M_PI*ypar[1]*ycos[1]*ycos[3]/L + M_PI*ypar[2]*ycos[1]*ycos[5]/L + 3*M_PI*ypar[9]*ycos[3]*ycos[7]/L + 3*M_PI*ycos[3]*ycos[7]*ycos[5]*ypar[11]/L + M_PI*ypar[3]*ycos[1]*ycos[3]*ycos[5]/L + 2*M_PI*ypar[4]*ycos[0]*ycos[2]*ycos[4]/L + 2*M_PI*ypar[5]*ycos[2]*ycos[6]*ycos[4]/L + 2*M_PI*ypar[6]*ycos[0]*ycos[8]*ycos[4]/L + 2*M_PI*ypar[7]*ycos[6]*ycos[8]*ycos[4]/L;
   return Pi2_2_result;
}

double Pi2_3(double *ypar,double *ycos, double *ysin, double L) {
   double Pi2_3_result;
   Pi2_3_result = -2*M_PI*ysin[1]*ysin[5]*ypar[2]/L - 2*M_PI*ysin[5]*ysin[7]*ypar[10]/L - 2*M_PI*ysin[1]*ysin[5]*ypar[3]*ycos[3]/L - M_PI*ysin[2]*ysin[4]*ypar[4]*ycos[0]/L - M_PI*ysin[2]*ysin[4]*ypar[5]*ycos[6]/L - 3*M_PI*ysin[4]*ysin[8]*ypar[6]*ycos[0]/L - 3*M_PI*ysin[4]*ysin[8]*ypar[7]*ycos[6]/L - 2*M_PI*ysin[5]*ysin[7]*ycos[3]*ypar[11]/L;
   return Pi2_3_result;
}

double Pi3_1(double *ypar,double *ycos, double *ysin, double L) {
   double Pi3_1_result;
   Pi3_1_result = -2*M_PI*ysin[2]*ysin[3]*ypar[2]/L - 2*M_PI*ysin[3]*ysin[8]*ypar[10]/L - M_PI*ysin[0]*ysin[5]*ypar[4]*ycos[1]/L - M_PI*ysin[0]*ysin[5]*ypar[5]*ycos[7]/L - 2*M_PI*ysin[2]*ysin[3]*ypar[3]*ycos[4]/L - 2*M_PI*ysin[3]*ysin[8]*ycos[4]*ypar[11]/L - 3*M_PI*ysin[6]*ysin[5]*ypar[6]*ycos[1]/L - 3*M_PI*ysin[6]*ysin[5]*ypar[7]*ycos[7]/L;
   return Pi3_1_result;
}

double Pi3_2(double *ypar,double *ycos, double *ysin, double L) {
   double Pi3_2_result;
   Pi3_2_result = -2*M_PI*ysin[2]*ysin[4]*ypar[1]/L - 2*M_PI*ysin[4]*ysin[8]*ypar[9]/L - M_PI*ysin[1]*ysin[5]*ypar[4]*ycos[0]/L - M_PI*ysin[1]*ysin[5]*ypar[6]*ycos[6]/L - 2*M_PI*ysin[2]*ysin[4]*ypar[3]*ycos[3]/L - 2*M_PI*ysin[4]*ysin[8]*ycos[3]*ypar[11]/L - 3*M_PI*ysin[5]*ysin[7]*ypar[5]*ycos[0]/L - 3*M_PI*ysin[5]*ysin[7]*ypar[7]*ycos[6]/L;
   return Pi3_2_result;
}

double Pi3_3(double *ypar, double *ycos, double *ysin, double L) {
   double Pi3_3_result;
   Pi3_3_result = M_PI*ypar[0]*ycos[2]/L + 3*M_PI*ypar[8]*ycos[8]/L + 3*M_PI*ycos[3]*ycos[8]*ypar[10]/L + M_PI*ypar[1]*ycos[2]*ycos[4]/L + M_PI*ypar[2]*ycos[3]*ycos[2]/L + 3*M_PI*ypar[9]*ycos[8]*ycos[4]/L + 3*M_PI*ycos[3]*ycos[8]*ycos[4]*ypar[11]/L + M_PI*ypar[3]*ycos[3]*ycos[2]*ycos[4]/L + 2*M_PI*ypar[4]*ycos[1]*ycos[0]*ycos[5]/L + 2*M_PI*ypar[5]*ycos[0]*ycos[7]*ycos[5]/L + 2*M_PI*ypar[6]*ycos[1]*ycos[6]*ycos[5]/L + 2*M_PI*ypar[7]*ycos[7]*ycos[6]*ycos[5]/L;
   return Pi3_3_result;
}

