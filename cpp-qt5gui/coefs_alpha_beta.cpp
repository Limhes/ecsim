
#include "coefs_alpha_beta.h"
//Alpha_N_M() and Beta_N_M() can be used to call all functions for forms (N,M)

double Alpha_N_M(int N, int M,int index, double gamma)
{
    switch (M)
         {
        case 1:return Alpha_N_1(N,index,gamma);
         case 2:return Alpha_N_2(N,index,gamma);
         case 3:return Alpha_N_3(N,index,gamma);
         case 4:return Alpha_N_4(N,index,gamma);
         case 5:return Alpha_N_5(N,index,gamma);
         default: return 0.;


         }
}
double Beta_N_M(int N, int M,int index, double gamma)
{
    switch (M)
         {
         case 1:return Beta_N_1(N,index,gamma);
         case 2:return Beta_N_2(N,index,gamma);
         case 3:return Beta_N_3(N,index,gamma);
         case 4:return Beta_N_4(N,index,gamma);
         case 5:return Beta_N_5(N,index,gamma);
         default: return 0.;


         }
}
double Alpha_N_2(int N,int index, double gamma)
{
    switch(N)
    {
    case 3:return Alpha_3_2(index,gamma);
    case 4:return Alpha_4_2(index,gamma);
    case 5:return Alpha_5_2(index,gamma);
    case 6:return Alpha_6_2(index,gamma);
    case 7:return Alpha_7_2(index,gamma);
    case 8:return Alpha_8_2(index,gamma);
    default: return 0.;
    }
}
double Alpha_N_1(int N,int index, double gamma)
{
    switch(N)
    {
    case 3:return Alpha_3_1(index,gamma);
    case 4:return Alpha_4_1(index,gamma);
    case 5:return Alpha_5_1(index,gamma);
   // case 6:return Alpha_6_1(index,gamma);
   // case 7:return Alpha_7_1(index,gamma);
   // case 8:return Alpha_8_1(index,gamma);
    default: return 0.;
    }
}

double Alpha_N_3(int N,int index, double gamma)
{
    switch(N)
    {
    case 3:return Alpha_3_3(index,gamma);
    case 4:return Alpha_4_3(index,gamma);
    case 5:return Alpha_5_3(index,gamma);
    case 6:return Alpha_6_3(index,gamma);

    default: return 0.;
    }
}
double Alpha_N_4(int N,int index, double gamma)
{
    switch(N)
    {

    case 4:return Alpha_4_4(index,gamma);
    case 5:return Alpha_5_4(index,gamma);
    case 6:return Alpha_6_4(index,gamma);

    default: return 0.;
    }
}

double Alpha_N_5(int N,int index, double gamma)
{
    switch(N)
    {

    case 5:return Alpha_5_5(index,gamma);
    case 6:return Alpha_6_5(index,gamma);

    default: return 0.;
    }
}
double Beta_N_2(int N,int index, double gamma)
{
    switch(N)
    {
    case 3:return Beta_3_2(index,gamma);
    case 4:return Beta_4_2(index,gamma);
    case 5:return Beta_5_2(index,gamma);
    case 6:return Beta_6_2(index,gamma);
    case 7:return Beta_7_2(index,gamma);

    default: return 0.;
    }
}
double Beta_N_3(int N,int index, double gamma)
{
    switch(N)
    {
    case 3:return Beta_3_3(index,gamma);
    case 4:return Beta_4_3(index,gamma);
    case 5:return Beta_5_3(index,gamma);
    case 6:return Beta_6_3(index,gamma);

    default: return 0.;
    }
}
double Beta_N_4(int N,int index, double gamma)
{
    switch(N)
    {

    case 4:return Beta_4_4(index,gamma);
    case 5:return Beta_5_4(index,gamma);
    case 6:return Beta_6_4(index,gamma);

    default: return 0.;
    }
}
double Beta_N_5(int N,int index, double gamma)
{
    switch(N)
    {


    case 5:return Beta_5_5(index,gamma);
    case 6:return Beta_6_5(index,gamma);

    default: return 0.;
    }
}
double Beta_N_1(int N,int index, double gamma)
{
    switch(N)
    {
    case 2:return (index==0)?-1.0:1.0; // added by Rene for 2-point current approximation
    case 3:return Beta_3_1(index,gamma);
    case 4:return Beta_4_1(index,gamma);
    case 5:return Beta_5_1(index,gamma);
    case 6:return Beta_6_1(index,gamma);


    default: return 0.;
    }
}

double Alpha_5_2(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        switch (index)
                {
                case -1:
                        result=2*gamma4*(3+4*gamma+3*gamma2+gamma3)/(1+gamma)/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                case 0:
                        result=2*(3-gamma-3*gamma2-3*gamma3-gamma4)/(1+gamma)/(1+gamma+gamma2);
                        break;
                case 1:
                        result=-2*(2+gamma-gamma2-2*gamma3-gamma4)/gamma2/(1+gamma)/(1+gamma);
                        break;
                case 2:
                        result=2*(2-gamma3)/gamma3/(1+gamma)/(1+gamma+gamma2);
                        break;
                case 3:
                        result=-2*(2-gamma2)/gamma3/(1+gamma)/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }
double Alpha_5_3(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        switch (index)
                {
                case -2:
                        result=-2*gamma4*gamma3*(2-gamma2)/(1+gamma)/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                case -1:
                        result=2*gamma3*(1+gamma)*(2-gamma)/(1+gamma+gamma2);
                        break;
                case 0:
                        result=2*(1-gamma-5*gamma2-gamma3+gamma4)/(1+gamma)/(1+gamma);
                        break;
                case 1:
                        result=-2*(1+gamma)*(1-2*gamma)/gamma/(1+gamma+gamma2);
                        break;
                case 2:
                        result=2*(1-2*gamma2)/gamma/(1+gamma)/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }
double Alpha_5_1(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4, gamma5, gamma6;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        gamma5=gamma4*gamma;
        gamma6=gamma5*gamma;
        switch (index)
                {
                case 0:
                        result=2*(6+9*gamma+9*gamma2+7*gamma3+3*gamma4+gamma5)/(1+gamma)/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                case 1:
                        result=-2*(3+3*gamma+4*gamma2+2*gamma3+gamma4)/gamma3/(1+gamma+gamma2);
                        break;
                case 2:
                        result=2*(3+4*gamma+5*gamma2+4*gamma3+2*gamma4+gamma5)/gamma5/(1+gamma)/(1+gamma);
                        break;
                case 3:
                        result=-2*(3+gamma+2*gamma2+gamma3)/gamma6/(1+gamma+gamma2);
                        break;
                case 4:
                        result=2*(3+4*gamma+3*gamma2+gamma3)/gamma6/(1+gamma)/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }
double Beta_5_1(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4, gamma5, gamma6;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        gamma5=gamma4*gamma;
        gamma6=gamma5*gamma;
        switch (index)
                {
                case 0:
                        result=-(4+5*gamma+7*gamma2+5*gamma3+3*gamma4+gamma5)/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                case 1:
                        result=(1+gamma)*(1+gamma2)/gamma3;
                        break;
                case 2:
                        result=-(1+gamma2)*(1+gamma+gamma2)/gamma5/(1+gamma);
                        break;
                case 3:
                        result=(1+gamma)*(1+gamma2)/gamma6/(1+gamma+gamma2);
                        break;
                case 4:
                        result=-1/gamma6/(1+gamma)/(1+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }
double Beta_5_2(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        switch (index)
                {
                case -1:
                        result=-gamma4/(1+gamma)/(1+gamma2);
                        break;
                case 0:
                        result=-(3+3*gamma+gamma2-gamma3-gamma4)/(1+gamma)/(1+gamma+gamma2);
                        break;
                case 1:
                        result=(1+gamma+gamma2)/gamma2/(1+gamma);
                        break;
                case 2:
                        result=-1/gamma3/(1+gamma);
                        break;
                case 3:
                        result=1/gamma3/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }
double Beta_5_3(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        switch (index)
                {
                case -2:
                        result=gamma3*gamma4/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                case -1:
                        result=-gamma3*(1+gamma)/(1+gamma+gamma2);
                        break;
                case 0:
                        result=-2*(1-gamma);
                        break;
                case 1:
                        result=(1+gamma)/gamma/(1+gamma+gamma2);
                        break;
                case 2:
                        result=-1/gamma/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }

double Beta_5_4(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        switch (index)
                {
                case -3:
                        result=-gamma4*gamma4*gamma/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                case -2:
                        result=gamma*gamma4/(1+gamma);
                        break;
                case -1:
                        result=-gamma2*(1+gamma+gamma2)/(1+gamma);
                        break;
                case 0:
                        result=-(1+gamma-gamma2-3*gamma3-3*gamma4)/(1+gamma)/(1+gamma+gamma2);
                        break;
                case 1:
                        result=1/(1+gamma+gamma2+gamma3);
                        break;
                default:
                        result=0;
                }
        return result;
        }
double Alpha_5_4(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        switch (index)
                {
                case -3:
                        result=2*gamma4*gamma4*gamma*(1-2*gamma2)/(1+gamma)/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                case -2:
                        result=-2*gamma*gamma4*(1-2*gamma3)/(1+gamma)/(1+gamma+gamma2);
                        break;
                case -1:
                        result=2*gamma2*(1+2*gamma+gamma2-gamma3-2*gamma4)/(1+gamma)/(1+gamma);
                        break;
                case 0:
                        result=-2*gamma*(1+3*gamma+3*gamma2+gamma3-3*gamma4)/(1+gamma)/(1+gamma+gamma2);
                        break;
                case 1:
                        result=2*gamma*(1+3*gamma+4*gamma2+3*gamma3)/(1+gamma)/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }
double Beta_5_5(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        switch (index)
                {
                case -4:
                        result=gamma4*gamma4*gamma2/(1+gamma)/(1+gamma2);
                        break;
                case -3:
                        result=-gamma3*gamma3*(1+gamma+gamma2+gamma3)/(1+gamma+gamma2);
                        break;
                case -2:
                        result=gamma3*(1+gamma2)*(1+gamma+gamma2)/(1+gamma);
                        break;
                case -1:
                        result=-gamma*(1+gamma)*(1+gamma2);
                        break;
                case 0:
                        result=gamma*(1+3*gamma+5*gamma2+7*gamma3+5*gamma4+4*gamma4*gamma)/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }
double Alpha_5_5(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        switch (index)
                {
                case -4:
                        result=2*gamma4*gamma4*gamma3*(1+3*gamma+4*gamma2+3*gamma3)/(1+gamma)/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                case -3:
                        result=-2*gamma4*gamma3*(1+2*gamma+gamma2+3*gamma3)/(1+gamma+gamma2);
                        break;
                case -2:
                        result=2*gamma4*(1+2*gamma+4*gamma2+5*gamma3+4*gamma4+3*gamma4*gamma)/(1+gamma)/(1+gamma);
                        break;
                case -1:
                        result=-2*gamma3*(1+2*gamma+4*gamma2+3*gamma3+3*gamma4)/(1+gamma+gamma2);
                        break;
                case 0:
                        result=2*gamma3*(1+3*gamma+7*gamma2+9*gamma3+9*gamma4+6*gamma4*gamma)/(1+gamma)/(1+gamma)/(1+gamma2)/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }
double Beta_3_2(int index, double gamma)
        {
        double result;
        double gamma2;
        gamma2=gamma*gamma;
        switch (index)
                {
                case -1:
                        result=-gamma2/(1+gamma);
                        break;
                case 0:
                        result=(gamma-1);
                        break;
                case 1:
                        result=1/(1+gamma);
                        break;
                default:
                        result=0;
                }
        return result;
        }
double Alpha_3_2(int index, double gamma)
        {
         double result;
        double gamma2;
        gamma2=gamma*gamma;
        switch (index)
                {
                case -1:
                        result=2*gamma2/(1+gamma);
                        break;
                case 0:
                        result=-2*gamma;
                        break;
                case 1:
                        result=2*gamma/(1+gamma);
                        break;
                default:
                        result=0;
                }
        return result;
       }
double Beta_3_1(int index, double gamma)
        {
        double result;
        switch (index)
                {
                case 0:
                        result=-(2+gamma)/(1+gamma);
                        break;
                case 1:
                        result=(1+gamma)/gamma;
                        break;
                case 2:
                        result=-1./gamma/(1+gamma);
                        break;
                default:
                        result=0;
                }
        return result;
        }
double Alpha_3_1(int index, double gamma)
        {
         double result;
         switch (index)
                {
                case 0:
                        result=2./(1+gamma);
                        break;
                case 1:
                        result=-2/gamma;
                        break;
                case 2:
                        result=2./gamma/(1+gamma);
                        break;
                default:
                        result=0;
                }
        return result;
       }
double Beta_3_3(int index, double gamma)
        {
        double result;
        double gamma2,gamma3;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        switch (index)
                {
                case -2:
                        result=gamma3/(1+gamma);
                        break;
                case -1:
                        result=-gamma*(1+gamma);
                        break;
                case 0:
                        result=gamma*(1+2*gamma)/(1+gamma);
                        break;
                default:
                        result=0;
                }
        return result;
        }
double Alpha_3_3(int index, double gamma)
        {
        double result;
        double gamma2,gamma3,gamma4;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
         switch (index)
                {
                case -2:
                        result=2*gamma4/(1+gamma);
                        break;
                case -1:
                        result=-2*gamma3;
                        break;
                case 0:
                        result=2*gamma3/(1+gamma);
                        break;
                default:
                        result=0;
                }
        return result;
       }

double Beta_4_1(int index, double gamma)
        {
        double result;
        double gamma2, gamma3;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;

        switch (index)
                {
                case 0:
                        result=-(3+4*gamma+3*gamma2+gamma3)/(1+gamma)/(1+gamma+gamma2);
                        break;
                case 1:
                        result=(1+gamma+gamma2)/gamma2;
                        break;
                case 2:
                        result=-(1+gamma+gamma2)/gamma3/(1+gamma);
                        break;
                case 3:
                        result=1/gamma3/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }

double Alpha_4_1(int index, double gamma)
        {
        double result;
        double gamma2, gamma3;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;

        switch (index)
                {
                case 0:
                        result=2*(3+2*gamma+gamma2)/(1+gamma)/(1+gamma+gamma2);
                        break;
                case 1:
                        result=-2*(2+2*gamma+gamma2)/gamma2/(1+gamma);
                        break;
                case 2:
                        result=2*(2+gamma+gamma2)/gamma3/(1+gamma);
                        break;
                case 3:
                        result=-2*(2+gamma)/gamma3/(1+gamma)/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }

double Beta_4_2(int index, double gamma)
        {
        double result;
        double gamma2, gamma3;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;

        switch (index)
                {
                case -1:
                        result=-gamma3/(1+gamma+gamma2);
                        break;
                case 0:
                        result=-(2-gamma2)/(1+gamma);
                        break;
                case 1:
                        result=1/gamma;
                        break;
                case 2:
                        result=-1/gamma/(1+gamma)/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }
double Alpha_4_2(int index, double gamma)
        {
        double result;
        double gamma2, gamma3;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;

        switch (index)
                {
               case -1:
                        result=2*gamma3*(2+gamma)/(1+gamma)/(1+gamma+gamma2);
                        break;
                case 0:
                        result=2*(1-2*gamma-gamma2)/(1+gamma);
                        break;
                case 1:
                        result=-2*(1-gamma-gamma2)/gamma/(1+gamma);
                        break;
                case 2:
                        result=2*(1-gamma)/gamma/(1+gamma)/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                 }
        return result;
        }
double Beta_4_3(int index, double gamma)
        {
        double result;
        double gamma2, gamma3;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;

        switch (index)
                {
                case -2:
                        result=gamma2*gamma3/(1+gamma)/(1+gamma+gamma2);
                        break;
                case -1:
                        result=-gamma2;
                        break;
                case 0:
                        result=-(1-2*gamma2)/(1+gamma);
                        break;
                case 1:
                        result=1/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }

double Alpha_4_3(int index, double gamma)
        {
        double result;
        double gamma2, gamma3;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;

        switch (index)
                {
                case -2:
                        result=-2*gamma2*gamma3*(1-gamma)/(1+gamma)/(1+gamma+gamma2);
                        break;
                case -1:
                        result=2*gamma2*(1+gamma-gamma2)/(1+gamma);
                        break;
                case 0:
                        result=-2*gamma*(1+2*gamma-gamma2)/(1+gamma);
                        break;
                case 1:
                        result=2*gamma*(1+2*gamma)/(1+gamma)/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }

double Beta_4_4(int index, double gamma)
        {
        double result;
        double gamma2, gamma3;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;

        switch (index)
                {
                case -3:
                        result=-gamma3*gamma3/(1+gamma+gamma2);
                        break;
                case -2:
                        result=gamma3*(1+gamma+gamma2)/(1+gamma);
                        break;
                case -1:
                        result=-gamma*(1+gamma+gamma2);
                        break;
                case 0:
                        result=gamma*(1+3*gamma+4*gamma2+3*gamma3)/(1+gamma)/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }

double Alpha_4_4(int index, double gamma)
        {
        double result;
        double gamma2, gamma3;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;

         switch (index)
                {
                case -3:
                        result=-2*gamma3*gamma3*gamma*(1+2*gamma)/(1+gamma)/(1+gamma+gamma2);
                        break;
                case -2:
                        result=2*gamma*gamma3*(1+gamma+2*gamma2)/(1+gamma);
                        break;
                case -1:
                        result=-2*gamma3*(1+2*gamma+2*gamma2)/(1+gamma);
                        break;
                case 0:
                        result=2*gamma3*(1+2*gamma+3*gamma2)/(1+gamma)/(1+gamma+gamma2);
                        break;
                default:
                        result=0;
                }
        return result;
        }


double Alpha_6_2(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4, gamma5;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        gamma5=gamma4*gamma;
        switch (index)
                {
                case -1:
                        result= 2*(gamma5+3*gamma4+5*gamma3+7*gamma2+5*gamma+4)*gamma5
                                  /(1+gamma)/(1+gamma+gamma2)/(gamma2+gamma3+gamma4+gamma+1)/(gamma2+1);
                        break;
                case 0:
                        result= -2*(gamma2*gamma5+4*gamma*gamma5+7*gamma5+9*gamma4+5*gamma3-5*gamma-6)
                                   /(gamma2+1)/(1+gamma+gamma2)/(1+gamma)/(1+gamma);
                        break;
                case 1:
                        result=2*(gamma*gamma5+2*gamma5+2*gamma4+gamma3-2*gamma2-2*gamma-3)
                                  /(1+gamma)/(1+gamma+gamma2)/gamma3;
                        break;
                case 2:
                        result=-2*(gamma*gamma5+gamma5+gamma4-gamma3-3*gamma2-3*gamma-3)
                                   /(1+gamma+gamma2)/(1+gamma)/(1+gamma)/gamma5;
                        break;
                case 3:
                        result=2*(gamma4-gamma2-3)
                                  /(1+gamma+gamma2)/(1+gamma)/(1+gamma2)/gamma/gamma5;
                        break;

                case 4: result= -2*(gamma4+gamma3-gamma2-3*gamma-3)
                                    /(1+gamma+gamma2)/(gamma2+gamma3+gamma4+gamma+1)/(1+gamma2)/(1+gamma)/(1+gamma)/gamma/gamma5;
                        break;
                default:
                        result=0;
                }
        return result;
        }

double Alpha_7_2(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4, gamma5, gamma6, gamma7, gamma8,gamma9;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        gamma5=gamma4*gamma;
        gamma6=gamma5*gamma;
        gamma7=gamma6*gamma;
        gamma8=gamma7*gamma;
        gamma9=gamma8*gamma;
        switch (index)
                {
                case -1:
                        result= 2*gamma6*(gamma9+4*gamma8+9*gamma7+16*gamma6+22*gamma5+26*gamma4+24*gamma3+19*gamma2+11*gamma+5)
                                   /(gamma2-gamma+1)/(1+gamma+gamma2+gamma3+gamma4)/(gamma2+1)/(1+gamma)/(1+gamma)/(1+gamma+gamma2)/(1+gamma+gamma2);
                        break;
                case 0:
                        result= -2*(gamma2*gamma9+5*gamma*gamma9+12*gamma9+21*gamma8+27*gamma7+27*gamma6+17*gamma5+gamma4-13*gamma3-20*gamma2-19*gamma-10)
                                   /(1+gamma+gamma2)/(1+gamma+gamma2+gamma3+gamma4)/(gamma2+1)/(1+gamma)/(1+gamma);
                        break;
                case 1:
                        result=2*(gamma*gamma9+3*gamma9+5*gamma8+6*gamma7+4*gamma6-6*gamma4-9*gamma3-10*gamma2-7*gamma-4)
                                  /(gamma2+1)/(1+gamma+gamma2)/(1+gamma)/(1+gamma)/gamma4;
                        break;
                case 2:
                        result= -2*(gamma*gamma9+2*gamma9+3*gamma8+2*gamma7-gamma6-6*gamma5-11*gamma4-13*gamma3-12*gamma2-8*gamma-4)
                                   /(1+gamma)/(1+gamma)/(1+gamma+gamma2)/(1+gamma+gamma2)/gamma7;
                        break;
                case 3:
                        result= 2*(gamma8+gamma7-3*gamma4-5*gamma3-5*gamma2-4*gamma-4)
                                   /(gamma2+1)/(1+gamma+gamma2)/(1+gamma)/(1+gamma)/gamma9;
                        break;

                case 4: result= -2*(gamma8+2*gamma7+gamma6-2*gamma5-5*gamma4-7*gamma3-9*gamma2-8*gamma-4)
                                   /(1+gamma+gamma2)/(1+gamma+gamma2+gamma3+gamma4)/(1+gamma2)/(1+gamma)/(1+gamma)/gamma/gamma9;
                        break;

                case 5: result= 2*(gamma6+gamma5-2*gamma3-5*gamma2-4*gamma-4)
                                   /(gamma2-gamma+1)/(1+gamma+gamma2+gamma3+gamma4)/(gamma2+1)/(1+gamma)/(1+gamma)/(1+gamma+gamma2)/(1+gamma+gamma2)/gamma/gamma9;
                        break;
                default:
                        result=0;
                }
        return result;
        }
double Beta_7_2(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4, gamma5, gamma6, gamma7, gamma8,gamma9,gamma10;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        gamma5=gamma4*gamma;
        gamma6=gamma5*gamma;
        gamma7=gamma6*gamma;
        gamma8=gamma7*gamma;
        gamma9=gamma8*gamma;
        gamma10=gamma9*gamma;
        switch (index)
                {
                case -1:
                        result= -gamma6/(gamma3+1+gamma4+gamma+gamma5+gamma2);
                        break;
                case 0:
                        result=(gamma10+2*gamma9+2*gamma8-5*gamma6-11*gamma5-17*gamma4-18*gamma3-16*gamma2-10*gamma-5)/
                            (gamma7+3*gamma6+5*gamma5+6*gamma4+6*gamma3+5*gamma2+3*gamma+1)/(gamma2+1);
                        break;
                case 1:
                        result=(1+gamma+gamma2+gamma3+gamma4)/(gamma+1)/gamma4;
                        break;
                case 2:
                    result= -(gamma6+gamma5+2*gamma4+2*gamma3+2*gamma2+gamma+1)/
                            gamma7/(2*gamma+2*gamma2+gamma3+1);
                        break;
                case 3:
                        result= (1+gamma+gamma2+gamma3+gamma4)/
                            gamma9/(2*gamma+2*gamma2+gamma3+1);
                        break;

                case 4: result= -1/gamma10/(1+gamma+gamma2+gamma3);
                        break;

                case 5: result= 1/gamma10/(2*gamma+3*gamma2+5*gamma5+4*gamma3+5*gamma4+4*gamma6+3*gamma7+2*gamma8+gamma9+1);
                        break;
                default:
                        result=0;
                }
        return result;
        }



double Alpha_8_2(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4, gamma5, gamma6, gamma7, gamma8,gamma9, gamma10,gamma11,gamma12,gamma13,gamma14,gamma15;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        gamma5=gamma4*gamma;
        gamma6=gamma5*gamma;
        gamma7=gamma6*gamma;
        gamma8=gamma7*gamma;
        gamma9=gamma8*gamma;
        gamma10=gamma9*gamma;
        gamma11=gamma10*gamma;
        gamma12=gamma11*gamma;
        gamma13=gamma12*gamma;
        gamma14=gamma13*gamma;
        gamma15=gamma14*gamma;
        switch (index)
                {
                case -1:
                        result= 2*(gamma11+3*gamma10+6*gamma9+11*gamma8+15*gamma7+21*gamma6+21*gamma5+23*gamma4+18*gamma3+15*gamma2+7*gamma+6)*gamma7
                                   /(1+gamma)/(gamma2+1)/(1+gamma+gamma2)/(gamma2-gamma+1)/(1+gamma+gamma2+gamma3+gamma4)/(gamma2+gamma3+gamma4+gamma5+gamma6+gamma+1);
                        break;
                case 0:
                        result= -2*(gamma15+5*gamma14+13*gamma13+26*gamma12+40*gamma11+54*gamma10+58*gamma9+51*gamma8+31*gamma7+3*gamma6-26*gamma5-46*gamma4-50*gamma3-46*gamma2-29*gamma-15)
                                   /(gamma2-gamma+1)/(1+gamma+gamma2+gamma3+gamma4)/(gamma2+1)/(1+gamma)/(1+gamma)/(1+gamma+gamma2)/(1+gamma+gamma2);
                        break;
                case 1:
                        result=  2*(gamma12+2*gamma11+3*gamma10+4*gamma9+3*gamma8+2*gamma7-3*gamma6-4*gamma5-9*gamma4-8*gamma3-9*gamma2-4*gamma-5)
                                    /(1+gamma)/(gamma2+1)/(1+gamma+gamma2+gamma3+gamma4)/gamma5;
                        break;
                case 2:
                        result= -2*(gamma12+gamma11+2*gamma10+gamma9-3*gamma7-7*gamma6-9*gamma5-12*gamma4-11*gamma3-10*gamma2-5*gamma-5)
                                    /(gamma2+1)/(1+gamma+gamma2)/(1+gamma)/(1+gamma)/gamma9;
                        break;
                case 3:
                        result= 2*(gamma12+gamma11+gamma10+gamma9-2*gamma8-4*gamma7-9*gamma6-10*gamma5-14*gamma4-11*gamma3-11*gamma2-5*gamma-5)
                                   /(gamma2+1)/(1+gamma)/(1+gamma+gamma2)/(1+gamma+gamma2)/gamma12;
                        break;

                case 4: result= -2*(gamma10+gamma9-gamma7-2*gamma6-4*gamma5-8*gamma4-7*gamma3-6*gamma2-5*gamma-5)
                                   /(gamma2+1)/(1+gamma+gamma3+gamma3+gamma4)/(1+gamma)/(1+gamma)/gamma14;
                        break;

                case 5: result= 2*(gamma8-gamma5-3*gamma4-gamma3-6*gamma2-5)
                                   /(1+gamma)/(1+gamma+gamma2)/(gamma2-gamma+1)/(1+gamma+gamma2+gamma3+gamma4)/(gamma2+1)/gamma15;
                        break;
                case 6: result= -2*(gamma10+2*gamma9+2*gamma8-5*gamma6-11*gamma5-17*gamma4-18*gamma3-16*gamma2-10*gamma-5)
                                   /(gamma2-gamma+1)/(1+gamma+gamma2+gamma3+gamma4)/(gamma2+gamma3+gamma4+gamma5+gamma6+gamma+1)/(gamma2+1)/(1+gamma+gamma2)/(1+gamma+gamma2)/(1+gamma)/(1+gamma)/gamma15;
                        break;

                default:
                        result=0;
                }
        return result;
        }
double Beta_6_2(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4, gamma5, gamma6, gamma7;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        gamma5=gamma4*gamma;
        gamma6=gamma5*gamma;
        gamma7=gamma6*gamma;
        switch (index)
                {
                case -1:
                        result= -gamma5/(gamma2+gamma3+gamma4+gamma+1);
                        break;
                case 0:
                        result= (gamma6+gamma5-2*gamma3-5*gamma2-4*gamma-4)/(1+gamma)/(gamma2+1)/(1+gamma+gamma2);
                        break;
                case 1:
                        result=(gamma2+1)/gamma3;
                        break;
                case 2:
                        result= -(gamma2+1)/(1+gamma)/gamma5;
                        break;
                case 3:
                        result=1/(1+gamma+gamma2)/gamma6;
                        break;

                case 4: result=-1/(gamma7+2*gamma6+3*gamma5+4*gamma4+4*gamma3+3*gamma2+2*gamma+1)/gamma6;
                        break;
                default:
                        result=0;
                }
        return result;
        }


double Beta_6_1(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4, gamma5, gamma6,gamma7,gamma8,gamma9,gamma10;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        gamma5=gamma4*gamma;
        gamma6=gamma5*gamma;
        gamma7=gamma6*gamma;
        gamma8=gamma7*gamma;
        gamma9=gamma8*gamma;
        gamma10=gamma9*gamma;

        switch (index)
                {
                case 0:
                        result=-(gamma9+4*gamma8+9*gamma7+16*gamma6+22*gamma5+26*gamma4+24*gamma3+19*gamma2+11*gamma+5)/
                            (1+gamma+gamma2+gamma3+gamma4)/(gamma3+gamma2+gamma+1)/(1+gamma+gamma2);
                        break;

                case 1:
                        result=(1+gamma+gamma2+gamma3+gamma4)/gamma4;
                        break;
                case 2:
                        result=-(gamma6+gamma5+2*gamma4+2*gamma3+2*gamma2+gamma+1)/
                            (1+gamma)/gamma7;
                        break;
                case 3:
                        result=(gamma6+gamma5+2*gamma4+2*gamma3+2*gamma2+gamma+1)/
                            (1+gamma+gamma2)/gamma9;
                        break;
                case 4:
                        result=-(1+gamma+gamma2+gamma3+gamma4)/(gamma3+gamma2+gamma+1)/gamma10;
                        break;
                case 5: result=1./(1+gamma+gamma2+gamma3+gamma4)/gamma10;
                        break;

                default:
                        result=0;
                }
        return result;
        }
double Beta_6_3(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4, gamma5, gamma6,gamma7,gamma8,gamma9;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        gamma5=gamma4*gamma;
        gamma6=gamma5*gamma;
        gamma7=gamma6*gamma;
        gamma8=gamma7*gamma;
        gamma9=gamma8*gamma;


        switch (index)
                {
                case -2:
                        result=gamma9/(3*gamma5+4*gamma4+gamma7+2*gamma6+1+2*gamma+3*gamma2+4*gamma3);
                        break;

                case -1:
                        result=-gamma4/(gamma2+1);
                        break;
                case 0:
                        result=(2*gamma3-3)/(1+gamma+gamma2);
                        break;
                case 1:
                        result=1/gamma2;
                        break;
                case 2:
                     result=-1/gamma3/(1+gamma)/(1+gamma2);
                        break;
                case 3: result=1/gamma3/(gamma8+2*gamma7+4*gamma6+5*gamma5+6*gamma4+5*gamma3+4*gamma2+2*gamma+1);
                        break;

                default:
                        result=0;
                }
        return result;
        }
double Alpha_6_3(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4, gamma5, gamma6,gamma7,gamma8,gamma9,gamma10;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        gamma5=gamma4*gamma;
        gamma6=gamma5*gamma;
        gamma7=gamma6*gamma;
        gamma8=gamma7*gamma;
        gamma9=gamma8*gamma;
        gamma10=gamma9*gamma;


        switch (index)
                {
                case -2:
                        result=2*(gamma4+gamma3-gamma2-3*gamma-3)*gamma9/
                            (20*gamma4+15*gamma3+9*gamma8+1+9*gamma2+4*gamma+gamma10+4*gamma9+15*gamma7+22*gamma5+20*gamma6);
                        break;

                case -1:
                    result=-2*(gamma3-gamma2-gamma-3)*gamma4/
                            (gamma4+gamma3+2*gamma2+gamma+1);
                        break;
                case 0:
                    result=2*(gamma6-5*gamma4-9*gamma3-7*gamma2+2*gamma+3)/
                            (1+gamma+gamma2)/(gamma2+2*gamma+1);
                        break;
                case 1:
                        result=2*(2*gamma3+gamma2+gamma-2)/
                            gamma2/(1+gamma+gamma2);
                        break;
                case 2:
                        result=-2*(2*gamma4+2*gamma3+gamma2-2*gamma-2)/
                            gamma3/(gamma6+6*gamma3+5*gamma2+3*gamma5+5*gamma4+1+3*gamma);
                        break;
                case 3: result=4*(gamma-1)/
                            gamma3/(gamma8+2*gamma7+4*gamma6+5*gamma5+6*gamma4+5*gamma3+4*gamma2+2*gamma+1);
                        break;

                default:
                        result=0;
                }
        return result;
        }


double Beta_6_4(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4, gamma5, gamma6,gamma7,gamma8,gamma9,gamma10,gamma12;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        gamma5=gamma4*gamma;
        gamma6=gamma5*gamma;
        gamma7=gamma6*gamma;
        gamma8=gamma7*gamma;
        gamma9=gamma8*gamma;
        gamma10=gamma9*gamma;
        gamma12=gamma10*gamma2;


        switch (index)
                {
                case -3:
                        result=-gamma12/(1+gamma+gamma2)/
                            (gamma6+gamma5+2*gamma4+2*gamma3+2*gamma2+gamma+1);
                        break;

                case -2:
                        result=gamma7/(1+gamma+gamma2+gamma3);
                        break;
                case -1:
                        result=-gamma3;
                        break;
                case 0:
                        result=(3*gamma3-2)/(1+gamma+gamma2);
                        break;
                case 1:
                        result=1/gamma/(gamma2+1);
                        break;
                case 2: result=-1/gamma/(3*gamma5+4*gamma4+gamma7+2*gamma6+1+2*gamma+3*gamma2+4*gamma3);
                        break;

                default:
                        result=0;
                }
        return result;
        }

double Alpha_6_4(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4, gamma5, gamma6,gamma7,gamma8,gamma9,gamma10,gamma12;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        gamma5=gamma4*gamma;
        gamma6=gamma5*gamma;
        gamma7=gamma6*gamma;
        gamma8=gamma7*gamma;
        gamma9=gamma8*gamma;
        gamma10=gamma9*gamma;
        gamma12=gamma10*gamma2;


        switch (index)
                {
                case -3:
                        result=-4*gamma12*(gamma-1)/
                            (gamma8+2*gamma7+4*gamma6+5*gamma5+6*gamma4+5*gamma3+4*gamma2+2*gamma+1);
                        break;

                case -2:
                        result=2*gamma7*(2*gamma4+2*gamma3-gamma2-2*gamma-2)/
                            (1+gamma)/(1+2*gamma+3*gamma2+3*gamma3+2*gamma4+gamma5);
                        break;
                case -1:
                        result=-2*(2*gamma3-gamma2-gamma-2)*gamma3/
                            (1+gamma+gamma2);
                        break;
                case 0:
                    result=2*(3*gamma6+2*gamma5-7*gamma4-9*gamma3-5*gamma2+1)/
                            (1+2*gamma+2*gamma2+gamma3)/(1+gamma);
                        break;
                case 1:
                        result=2*(3*gamma3+gamma2+gamma-1)/
                            gamma/(gamma4+gamma3+2*gamma2+gamma+1);
                        break;
                case 2: result=-2*(3*gamma4+3*gamma3+gamma2-gamma-1)/
                            gamma/(1+gamma)/(gamma9+9*gamma6+11*gamma5+3*gamma8+6*gamma7+3*gamma+6*gamma2+9*gamma3+11*gamma4+1);
                        break;

                default:
                        result=0;
                }
        return result;
        }

double Beta_6_5(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4, gamma5, gamma6,gamma7,gamma8,gamma9,gamma10,gamma12,gamma14;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        gamma5=gamma4*gamma;
        gamma6=gamma5*gamma;
        gamma7=gamma6*gamma;
        gamma8=gamma7*gamma;
        gamma9=gamma8*gamma;
        gamma10=gamma9*gamma;
        gamma12=gamma10*gamma2;
        gamma14=gamma12*gamma2;


        switch (index)
                {
                case -4:
                        result=gamma14/(1+gamma+gamma2+gamma3)/
                            (1+gamma+gamma2+gamma3+gamma4);
                        break;

                case -3:
                        result=-gamma9/(1+gamma+gamma2);
                        break;
                case -2:
                        result=(gamma2+1)*gamma5/(1+gamma);
                        break;
                case -1:
                        result=-(gamma2+1)*gamma2;
                        break;
                case 0:
                        result=(4*gamma6+4*gamma5+5*gamma4+2*gamma3-gamma-1)/
                                (1+gamma)/(1+gamma+gamma2)/(gamma2+1);
                        break;
                case 1: result=1/(1+gamma+gamma2+gamma3+gamma4);
                        break;

                default:
                        result=0;
                }
        return result;
        }

double Alpha_6_5(int index, double gamma)
        {
        double result;
        double gamma2, gamma3,gamma4, gamma5, gamma6,gamma7,gamma8,gamma9,gamma10,gamma12,gamma14;
        gamma2=gamma*gamma;
        gamma3=gamma2*gamma;
        gamma4=gamma3*gamma;
        gamma5=gamma4*gamma;
        gamma6=gamma5*gamma;
        gamma7=gamma6*gamma;
        gamma8=gamma7*gamma;
        gamma9=gamma8*gamma;
        gamma10=gamma9*gamma;
        gamma12=gamma10*gamma2;
        gamma14=gamma12*gamma2;


        switch (index)
                {
                case -4:
                        result=2*gamma14*(3*gamma4+3*gamma3+gamma2-gamma-1)/
                            (1+gamma+gamma2)/(5*gamma6+3*gamma7+gamma8+7*gamma5+8*gamma4+7*gamma3+5*gamma2+3*gamma+1);
                        break;

                case -3:
                        result=-2*(3*gamma4+gamma2-1)*gamma9/
                            (1+gamma+gamma2)/(1+gamma+gamma2+gamma3);
                        break;
                case -2:
                        result=2*(3*gamma6+3*gamma5+3*gamma4+gamma3-gamma2-gamma-1)*gamma5/
                            (1+2*gamma+gamma2)/(1+gamma+gamma2);
                        break;
                case -1:
                        result=-2*gamma2*(3*gamma6+2*gamma5+2*gamma4-gamma3-2*gamma2-2*gamma-1)/
                            (1+gamma+gamma2)/(1+gamma);
                        break;
                case 0:
                        result=2*gamma*(6*gamma7+5*gamma6-5*gamma4-9*gamma3-7*gamma2-4*gamma-1)/
                            (1+gamma+gamma2+gamma3)/(1+gamma+gamma2)/(1+gamma);
                        break;
                case 1: result=2*(4*gamma5+5*gamma4+7*gamma3+5*gamma2+3*gamma+1)*gamma/
                            (1+gamma+gamma2)/(3*gamma5+4*gamma4+gamma7+2*gamma6+1+2*gamma+3*gamma2+4*gamma3);
                        break;

                default:
                        result=0;
                }
        return result;
        }
