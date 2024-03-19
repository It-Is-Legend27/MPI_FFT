// Assignment #2
// Name: Angel Badillo and Chad Callender
// Parallel Programming
// Date: 03/29/24
//***************************************************************
//***************************************************************
// How to run:
// This program is to be run on the TACC cluster using the SBATCH
// shell script named "AngelBadilloA2Script".
// The command to be run in the bash terminal is:
// sbatch AngelBadilloA2Script
//
// Description:
// Text here.
//*****************************************************************
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define N_SAMPLES 1024
#define PI 3.141592653589793

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

//*******************************************************************
// CmplxAdd
// @param X The augend, a CmplxNum
// @param Y The addend, a CmplxNum
// @brief
// Calculates the addition of complex numbers (CmplxNum) X and Y, then
// returns sum as CmplxNum.
// @return The sum as struct CmplxNum
//********************************************************************
struct CmplxNum CmplxAdd(struct CmplxNum X, struct CmplxNum Y);

//*******************************************************************
// CmplxSub
// @param X The subtrahend, a CmplxNum
// @param Y The minuend, a CmplxNum
// @brief
// Calculates the subtraction of complex numbers (CmplxNum) X and Y, then
// returns the difference as CmplxNum.
// @return The difference as struct CmplxNum
//********************************************************************
struct CmplxNum CmplxSub(struct CmplxNum X, struct CmplxNum Y);

//*******************************************************************
// CmplxMult
// @param X The multiplicand, a CmplxNum
// @param Y The multiplier, a CmplxNum
// @brief
// Calculates the multiplication of complex numbers (CmplxNum) X and
// Y, then returns product as CmplxNum.
// @return The product as struct CmplxNum
//********************************************************************
struct CmplxNum CmplxMult(struct CmplxNum X, struct CmplxNum Y);

//*******************************************************************
// partial_FFT
// @param XR an array of doubles of size 8192 containing real parts
// of output for N no. of FTT coefficients.
// @param XI an array of doubles of size 8192 containing imaginary
// parts of output for N no. of FFT coefficients.
// @param R an array of doubles of size 8192 containing real part of
// samples for the function x(n).
// @param I an array of doubles of size 8192 containing imaginary
// part of samples for the function x(n)
// @brief
// Implements the Cooley-Turkey FFT algorithm (AKA Radix-2).
// Computes the output for a total of N no. FFT coefficients with N no. of
// samples."Returns" (or modifies via reference) array of doubles XR and XI to
// have the output for N no. of FFT coefficients X(0) to X(N-1).
// @return void
//********************************************************************
void partial_FFT(double *XR, double *XI, double *R, double *I, int N);

int main(void)
{
    int rank, comm_sz; // Process rank and number of processes
    int local_n; // Number of coefficients to be computed by each process

    double *R = NULL, *I = NULL;   // Input values for samples
    double *XR = NULL, *XI = NULL; // Output values for coefficients

    double *local_XR_l, *local_XI_l; // Locally computed coefficients [0, N/2 - 1]
    double *local_XR_r, *local_XI_r; // Locally compute coefficients [N/2, N-1]

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Number of coefficients to be computed by each process
    local_n = N_SAMPLES / comm_sz;

    // All proccesses require samples
    R = (double *)calloc(N_SAMPLES, sizeof(double));
    I = (double *)calloc(N_SAMPLES, sizeof(double));
    
    // Set values for first 8 samples, the rest will be 0's.
    memcpy(R, (double[]){3.6, 2.9, 5.6, 4.8, 3.3, 5.9, 5, 4.3}, sizeof(double)* 8);
    memcpy(I, (double[]){2.6, 6.3, 4, 9.1, 0.4, 4.8, 2.6, 4.1}, sizeof(double)* 8);
    
    // Arrays to store coefficients stored by each process
    local_XR_l = (double *)calloc(local_n / 2, sizeof(double));
    local_XI_l = (double *)calloc(local_n / 2, sizeof(double));
    local_XR_r = (double *)calloc(local_n / 2, sizeof(double));
    local_XI_r = (double *)calloc(local_n / 2, sizeof(double));

    // Process 0 does this
    if (!rank)
    {
        XR = (double *)calloc(N_SAMPLES, sizeof(double));
        XI = (double *)calloc(N_SAMPLES, sizeof(double));
    }



    MPI_Finalize();
    return 0;
}

struct CmplxNum CmplxAdd(struct CmplxNum X, struct CmplxNum Y)
{
    struct CmplxNum Z = {.a = X.a + Y.a, .bi = X.bi + Y.bi};
    return Z;
}

struct CmplxNum CmplxSub(struct CmplxNum X, struct CmplxNum Y)
{
    struct CmplxNum Z = {.a = X.a - Y.a, .bi = X.bi - Y.bi};
    return Z;
}

struct CmplxNum CmplxMult(struct CmplxNum X, struct CmplxNum Y)
{
    // X = a + bi
    // Y = c + di
    // XY = (a+bi)(c+di) = (ac−bd) + (ad+bc)i
    struct CmplxNum Z = {.a = (X.a * Y.a) - (X.bi * Y.bi), .bi = (X.a * Y.bi) + (X.bi * Y.a)};
    return Z;
}

void partial_FFT(double *XR, double *XI, double *R, double *I, int N)
{
    for (int k = 0; k < N / 2 - 1; k++)
    {
        // Calculate twiddle factor
        struct CmplxNum tFactor = {.a = cos(2 * PI * k / N), .bi = -sin(2 * PI * k / N)};

        // Initialize variables for calculation of even and odd parts for calculation of X(k) and X(k + N/2)
        struct CmplxNum evenPart = {.a = 0, .bi = 0};
        struct CmplxNum oddPart = {.a = 0, .bi = 0};

        // Calculate the sums of the even and odd parts for calculation of X(k) and X(k + N/2)
        for (int m = 0; m <= N / 2 - 1; m++)
        {
            // E(k)
            struct CmplxNum funcX2m = {.a = R[2 * m], .bi = I[2 * m]};
            struct CmplxNum eulerPart = {.a = cos(2 * PI * 2 * m * k / N), .bi = -sin(2 * PI * 2 * m * k / N)};
            struct CmplxNum resEven = CmplxMult(funcX2m, eulerPart);
            evenPart = CmplxAdd(evenPart, resEven);

            // O(k)
            struct CmplxNum funcX2mP1 = {.a = R[(2 * m) + 1], .bi = I[(2 * m) + 1]};
            struct CmplxNum resOdd = CmplxMult(funcX2mP1, eulerPart);
            oddPart = CmplxAdd(oddPart, resOdd);
        }

        // Adding even part for E(k)
        XR[k] = evenPart.a;
        XI[k] = evenPart.bi;

        // Adding even part for E(k + N/2)
        XR[k + N / 2] = evenPart.a;
        XI[k + N / 2] = evenPart.bi;

        // Adjusting odd parts by twiddle factor
        struct CmplxNum temp = CmplxMult(tFactor, oddPart);

        // Adding odd part for O(k)
        XR[k] += temp.a;
        XI[k] += temp.bi;

        // Subtracting odd part for O(k + N/2)
        XR[k + N / 2] -= temp.a;
        XI[k + N / 2] -= temp.bi;
    }
}