#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int N = 8192;
const int pi = 3.141592653589793;
/**
 * @brief Represents a Complex number.
 * Contains real and imaginary part.
 *
 */
struct CmplxNum
{
    double a;  // real part of complex number
    double bi; // imaginary part of complex number
};

/**
 * CmplxAdd
 * @brief
 * Calculates the sum of CmplxNum X and Y, then returns sum as CmplxNum.
 * @param X Augend, a CmplxNum
 * @param Y Addend, a CmplxNum
 * @return Sum as struct CmplxNum
 */
struct CmplxNum CmplxAdd(struct CmplxNum X, struct CmplxNum Y)
{
    struct CmplxNum Z = {.a = X.a + Y.a, .bi = X.bi + Y.bi};
    return Z;
}

/**
 * CmplxSub
 * @brief
 * Calculates the difference from CmplxNum X and Y, then returns difference as CmplxNum.
 * @param X Subtrahend, as a  CmplxNum
 * @param Y Minuend, as a CmplxNum
 * @return Difference, as a struct CmplxNum
 */
struct CmplxNum CmplxSub(struct CmplxNum X, struct CmplxNum Y)
{
    struct CmplxNum Z = {.a = X.a - Y.a, .bi = X.bi - Y.bi};
    return Z;
}

/**
 * CmplxMult
 * @brief
 * Multiplies two CmplxNum together, then returns product as a CmplxNum.
 * @param X Multiplicand, as a CmplxNum
 * @param Y Multiplier, as a CmplxNum
 * @return struct CmplxNum
 */
struct CmplxNum CmplxMult(struct CmplxNum X, struct CmplxNum Y)
{
    // X = a + bi
    // Y = c + di
    // XY = (a+bi)(c+di) = (ac−bd) + (ad+bc)i
    struct CmplxNum Z = {.a = (X.a * Y.a) - (X.bi * Y.bi), .bi = (X.a * Y.bi) + (X.bi * Y.a)};
    return Z;
}

struct CmplxNum evenPartAtm(double R[], double I[], int m, int N, int k)
{
    struct CmplxNum funcX2n = {.a = R[2 * m], .bi = I[2 * m]};
    struct CmplxNum eulerPart = {.a = cos(2 * pi * 2 * m * k / N), .bi = -sin(2 * pi * 2 * m * k / N)};
    return CmplxMult(funcX2n, eulerPart);
}

void evenPartOfK(double XR[], double XI[], double R[], double I[], int M, int k)
{
    for (int n = 0; n <= (M / 2) - 1; n++)
    {
        struct CmplxNum result = evenPartAtm(R, I, n, M, k);
        printf("Even n=%d: %f + %fi\n", n, result.a, result.bi);
        XR[k] += result.a;
        XI[k] += result.bi;
    }
}

struct CmplxNum oddPartAtm(double R[], double I[], int m, int N, int k)
{
    struct CmplxNum funcX2n = {.a = R[(2 * m) + 1], .bi = I[(2 * m) + 1]};
    struct CmplxNum eulerPart = {.a = cos(2 * pi * 2 * m * k / N), .bi = -sin(2 * pi * 2 * m * k / N)};
    return CmplxMult(funcX2n, eulerPart);
}

struct CmplxNum twiddleFactor(int N, int k)
{
    struct CmplxNum tFactor = {.a = cos(2 * pi * k / N), .bi = -sin(2 * pi * k / N)};
    return tFactor;
}

void oddPartOfK(double XR[], double XI[], double R[], double I[], int N, int k)
{
    for (int n = 0; n <= (N / 2) - 1; n++)
    {
        struct CmplxNum result = CmplxMult(twiddleFactor(N, k), oddPartAtm(R, I, n, N, k));
        printf("Odd n=%d: %f + %fi\n", n, result.a, result.bi);
        XR[k] += result.a;
        XI[k] += result.bi;
    }
}

int main()
{
    // const double pi = 2*acos(0.0);
    // printf("%.15f\n", pi);
    double R[8] = {3.6, 2.9, 5.6, 4.8, 3.3, 5.9, 5, 4.3};
    double I[8] = {2.6, 6.3, 4, 9.1, 0.4, 4.8, 2.6, 4.1};
    double XR[8] = {0};
    double XI[8] = {0};

    evenPartOfK(XR, XI, R, I, 8, 1);
    oddPartOfK(XR, XI, R, I, 8, 1);

    printf("XR[1]: %.6f        XI[1]: %.6fi\n", XR[1], XI[1]);

    // 3.141592653589793

    // double XR[K];
    // double XI[K];

    // printf("==========================================================================\n");
    // printf("TOTAL PROCESSED SAMPLES: %d\n", N);
    // printf("==========================================================================\n");
    // for (int i = 0; i < K; ++i)
    // {
    //     printf("XR[%d]: %.6f          XI[%d]: %.6fi\n", i, XR[i], i, XI[i]);
    //     printf("==========================================================================\n");
    // }
}