#include <stdio.h>
#include <math.h>

void Comp(double *S, double *E, double *I, double *R, double *D, double alphaI, double alphaE, double B, double u, double k, double p, double gamma, double N, double dt)
{
    double dS = (-1 * ((alphaI * *S * *I) / N + (alphaE * *S * *E) / N) + gamma * *R);
    double dE = (1 * ((alphaI * *S * *I) / N + (alphaE * *S * *E) / N) - (k + p) * *E);
    double dI = (k * *E - B * *I - u * *I);
    double dR = (B * *I + p * *E - gamma * *R);
    double dD = (u * *I);

    double delta = fabs(dS) + fabs(dE) + fabs(dI) + fabs(dR) + fabs(dD);
    if (delta < 0.01)
    {
        return;
    }
    *S += dS * dt;
    *E += dE * dt;
    *I += dI * dt;
    *R += dR * dt;
    *D += dD * dt;
}

int main()
{
    double B = 0.999, gamma = 0.001, u = 0.0188, alphaE = 0.999, alphaI = 0.999;
    double S = 2798170.0 - 99.0 - 24.0, E = 99.0, I = 0.0, R = 24.0, D = 0.0;
    
    int N = 2798170;
    double k = 0.042;
    double p = 0.952;
    int dt = 1;
    int T = 68;
    

    printf("%10s%10s%10s%10s%10s\n", "S", "E", "I", "R", "D");


    for (int i = 0; i < T / dt; ++i)
    {
        printf("%10.2f%10.2f%10.2f%10.2f%10.2f\n", S, E, I, R, D);
        Comp(&S, &E, &I, &R, &D, alphaI, alphaE, B, u, k, p, gamma, N, dt);
    }

    return 0;
}