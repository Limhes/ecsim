/*
 * this file is a close to one-on-one copy from the supporting information of Molina et al. (10.1016/j.cplett.2015.11.011)
 * changes:
 * - parameter list changed so that N and M are of type std::size_t & index is of type int & gamma is of type double
 * 
 * analytical expressions for the coefficients can be found in 10.1016/j.electacta.2008.08.039
 */

#ifndef COEFS_ALPHA_BETA
#define COEFS_ALPHA_BETA

#include <cmath>

double Alpha_N_M(std::size_t N, std::size_t M,int index, double gamma);
double Beta_N_M(std::size_t N, std::size_t M,int index, double gamma);

double Alpha_N_1(std::size_t N,int index, double gamma);
double Alpha_N_2(std::size_t N,int index, double gamma);
double Alpha_N_3(std::size_t N,int index, double gamma);
double Alpha_N_4(std::size_t N,int index, double gamma);
double Alpha_N_5(std::size_t N,int index, double gamma);

double Beta_N_1(std::size_t N,int index, double gamma);
double Beta_N_2(std::size_t N,int index, double gamma);
double Beta_N_3(std::size_t N,int index, double gamma);
double Beta_N_4(std::size_t N,int index, double gamma);
double Beta_N_5(std::size_t N,int index, double gamma);

double Alpha_5_2(int index, double gamma);
double Alpha_5_3(int index, double gamma);
double Alpha_5_1(int index, double gamma);
double Beta_5_2(int index, double gamma);
double Beta_5_3(int index, double gamma);
double Beta_5_1(int index, double gamma);

double Alpha_5_4(int index, double gamma);
double Alpha_5_5(int index, double gamma);
double Beta_5_4(int index, double gamma);
double Beta_5_5(int index, double gamma);

double Alpha_3_1(int index, double gamma);
double Alpha_3_2(int index, double gamma);
double Beta_3_1(int index, double gamma);
double Beta_3_2(int index, double gamma);
double Alpha_3_3(int index, double gamma);
double Beta_3_3(int index, double gamma);

double Alpha_4_1(int index, double gamma);
double Alpha_4_2(int index, double gamma);

double Beta_4_1(int index, double gamma);
double Beta_4_2(int index, double gamma);
double Alpha_4_3(int index, double gamma);
double Beta_4_3(int index, double gamma);
double Alpha_4_4(int index, double gamma);
double Beta_4_4(int index, double gamma);

double Alpha_6_2(int index, double gamma);
double Alpha_7_2(int index, double gamma);
double Alpha_8_2(int index, double gamma);

double Beta_6_2(int index, double gamma);
double Beta_7_2(int index, double gamma);

double Beta_6_1(int index, double gamma);
double Beta_7_2(int index, double gamma);

double Alpha_6_3(int index, double gamma);
double Alpha_6_4(int index, double gamma);
double Alpha_6_5(int index, double gamma);

double Beta_6_3(int index, double gamma);
double Beta_6_4(int index, double gamma);
double Beta_6_5(int index, double gamma);


#endif // COEFS_ALPHA_BETA
