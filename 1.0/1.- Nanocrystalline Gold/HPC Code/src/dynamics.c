#include <math.h>
#include <complex.h>
#include <fftw3-mpi.h>
#include <gsl/gsl_randist.h>

#include "binary.h"

static double ipow (double x, int p)
{
    /* Integer power of a number x^p */
    if (p == 1) return x;
    else return x * ipow (x, p-1);
}

void normalize (state *s, fftw_complex *field)
{
    /* Normalize a fourier space field after fft */

    int ij;
    double norm_scale = 1.0 /(s->N * s->N);

    for (int i = 0; i < s->local_n1; i++)
    {
        for (int j = 0; j < s->N; j++)
        {
            ij = i*s->N + j;
            field[ij] *= norm_scale;
        }
    }

    MPI_Barrier (MPI_COMM_WORLD);

    return;
}

static double F_mix (double c, state *s)
{
    /* Free energy of mixing */
    double epsilon = -4.0 + s->epsilon0*(s->sigma - s->sigma0);

    return s->omega*(
                c * log(2.0 * c) + (1 - c) * log(2.0 - 2.0 * c) +
                0.5 * epsilon * ipow(c - 0.5, 2)
            );
}

static double dF_mixdc (double c, state *s)
{
    /* Ideal chemical potential of mixing */
    return s->omega * log(c / (1.0 - c));
}

static void noise(state *s)
{
    /* Gaussian Noise in Fourier Space */
    fftw_complex c_scale, n_scale;
    c_scale = I*sqrt(s->Mc / (s->dx * s->dx * s->dt));
    n_scale = I*sqrt(s->Mn / (s->dx * s->dx * s->dt));

    int ij;

    for (int i = 0; i < s->local_n1; i++)
    {
        for (int j = 0; j < s->N; j++)
        {
            ij = i*s->N + j;
            s->fxin[ij] = sqrt(s->k2[ij]) * n_scale * gsl_ran_gaussian_ziggurat(s->rng, 1.0) +
                      I * sqrt(s->k2[ij]) * n_scale * gsl_ran_gaussian_ziggurat(s->rng, 1.0);

            s->fxic[ij] = sqrt(s->k2[ij]) * c_scale * gsl_ran_gaussian_ziggurat(s->rng, 1.0) +
                      I * sqrt(s->k2[ij]) * c_scale * gsl_ran_gaussian_ziggurat(s->rng, 1.0);
        }
    }

    MPI_Barrier (MPI_COMM_WORLD);

    return;
}

static void set_nonlinear (state *s)
{
    /* Compute the nonlinear part of the equation of motion */
    double third = 1.0/3.0;
    int ij;

    for (int i = 0; i < s->local_n0; i++)
    {
        for (int j = 0; j < s->N; j++)
        {
            ij = i*2*((s->N>>1) + 1) + j;

            s->nnl[ij]  = ipow(s->n[ij], 2)*( -0.5 * s->eta + third * s->n[ij]);
            s->nnl[ij] += F_mix (s->c[ij], s);
            s->nnl[ij] += -exp(-ipow(s->c[ij] - 1.0, 2)/(2.0*ipow(s->alphac, 2))) * s->Cn[ij];

            s->cnl[ij]  = (s->n[ij] + 1.0) *  dF_mixdc (s->c[ij], s);
            s->cnl[ij] += 0.5 * s->n[ij] * (s->c[ij] - 1.0) / ipow(s->alphac, 2) *
                          exp(-ipow(s->c[ij] - 1.0, 2)/(2.0*ipow(s->alphac, 2))) *
                          s->Cn[ij];
        }
    }

    MPI_Barrier (MPI_COMM_WORLD);

    return;
}

static void calccorr (state *s)
{
    /* Calculates the correlation term of the equation of motion */
    int ij;
    double norm = 1.0/(s->N * s->N);

    for (int i = 0; i < s->local_n1; i++)
    {
        for (int j = 0; j < s->N; j++)
        {
            ij = i * s->N + j;
            s->fn[ij] *= norm;
            s->fCn[ij] = s->C[ij]*s->fn[ij];
        }
    }

    fftw_mpi_execute_dft_c2r (s->ifft_plan, s->fCn, s->Cn);

    return;
}

static void propagate (state *s)
{
    /* Propagate the fields forward in fourier space */
    double norm = 1.0 / (s->N * s->N);
    int ij;

    for (int i = 0; i < s->local_n1; i++)
    {
        for (int j = 0; j < s->N; j++)
        {
            ij = i * s->N + j;

            /* Do some normalization quickly */
            s->fc[ij] *= norm;
            s->fnnl[ij] *= norm;
            s->fcnl[ij] *= norm;
            s->fxin[ij] *= norm;
            s->fxic[ij] *= norm;

            s->fn[ij]  = s->Pn[ij] * s->fn[ij] + s->Qn[ij] * s->fnnl[ij] +  
                         s->Ln[ij] * s->fxin[ij];

            s->fc[ij]  = s->Pc[ij] * s->fc[ij] + s->Qc[ij] * s->fcnl[ij] + 
                         s->Lc[ij] * s->fxic[ij];
        }
    }

    MPI_Barrier (MPI_COMM_WORLD);

    return;
}

void step(state *s)
{
    /* Fourier transform fields */
    fftw_mpi_execute_dft_r2c (s->fft_plan, s->n, s->fn);
    fftw_mpi_execute_dft_r2c (s->fft_plan, s->c, s->fc);

    /* Calculate (C * n) and the nonlinear term */
    calccorr (s);
    set_nonlinear (s);

    /* Fourier transform the nonlinear terms */
    fftw_mpi_execute_dft_r2c (s->fft_plan, s->nnl, s->fnnl);
    fftw_mpi_execute_dft_r2c (s->fft_plan, s->cnl, s->fcnl);

    /* Make some noise and propagate the fields in fourier space */
    noise (s);
    propagate (s);

    /* Inverse fourier transform the propagated fields */
    fftw_mpi_execute_dft_c2r (s->ifft_plan, s->fn, s->n);
    fftw_mpi_execute_dft_c2r (s->ifft_plan, s->fc, s->c);

    /* Update the time and step #   */
    s->t += s->dt;
    s->step += 1;

    MPI_Barrier (MPI_COMM_WORLD);

    return;
}
