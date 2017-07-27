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
    
    /* Local alias for fields (for brevity sake) */
    double *nnl = s->nnl->real;
    double *cnl = s->nnl->real;
    double *n = s->n->real;
    double *c = s->c->real;
    double *Cn = s->Cn->real;

    for (int i = 0; i < s->local_n0; i++)
    {
        for (int j = 0; j < s->N; j++)
        {
            ij = i*2*((s->N>>1) + 1) + j;

            nnl[ij]  = ipow(n[ij], 2) * ( -0.5 * s->eta + third * n[ij]);
            nnl[ij] += F_mix (c[ij], s);
            nnl[ij] += -exp(-ipow(c[ij] - 1.0, 2) / (2.0 * ipow(s->alphac, 2))) * Cn[ij];

            cnl[ij]  = (n[ij] + 1.0) *  dF_mixdc (c[ij], s);
            cnl[ij] += 0.5 * n[ij] * (c[ij] - 1.0) / ipow(s->alphac, 2) *
                       exp(-ipow(c[ij] - 1.0, 2) / (2.0 * ipow(s->alphac, 2))) * Cn[ij];
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

    /* Local alias for fourier representations */
    fftw_complex *n = s->n->fourier;
    fftw_complex *Cn = s->Cn->fourier; 

    for (int i = 0; i < s->local_n1; i++)
    {
        for (int j = 0; j < s->N; j++)
        {
            ij = i * s->N + j;
            n[ij] *= norm;
            Cn[ij] = s->C[ij] * n[ij];
        }
    }

    ifft (s->ifft_plan, s->Cn);

    return;
}


static void propagate (state *s)
{
    /* Propagate the fields forward in fourier space */
    double norm = 1.0 / (s->N * s->N);
    int ij;

    /* Local alias for fields */
    fftw_complex *n = s->n->fourier;
    fftw_complex *c = s->c->fourier;
    fftw_complex *nnl = s->nnl->fourier;
    fftw_complex *cnl = s->cnl->fourier;

    for (int i = 0; i < s->local_n1; i++)
    {
        for (int j = 0; j < s->N; j++)
        {
            ij = i * s->N + j;

            /* Do some normalization quickly */
            c[ij] *= norm;
            nnl[ij] *= norm;
            cnl[ij] *= norm;
            s->fxin[ij] *= norm;
            s->fxic[ij] *= norm;

            /* Propagate in the Fourier domain! */
            n[ij] = s->Pn[ij] * n[ij] + s->Qn[ij] * nnl[ij] +  s->Ln[ij] * s->fxin[ij];
            c[ij]  = s->Pc[ij] * c[ij] + s->Qc[ij] * cnl[ij] + s->Lc[ij] * s->fxic[ij];
        }
    }

    MPI_Barrier (MPI_COMM_WORLD);

    return;
}


void step(state *s)
{
    /* Step the state forward in time one time step */
    
    /* Fourier transform fields */
    fft (s->fft_plan, s->n);
    fft (s->fft_plan, s->c);

    /* Calculate (C * n) and the nonlinear term */
    calccorr (s);
    set_nonlinear (s);

    /* Fourier transform the nonlinear terms */
    fft (s->fft_plan, s->nnl);
    fft (s->fft_plan, s->cnl);

    /* Make some noise and propagate the fields in fourier space */
    noise (s);
    propagate (s);

    /* Inverse fourier transform the propagated fields */
    ifft (s->ifft_plan, s->n);
    ifft (s->ifft_plan, s->c);

    /* Update the time and step #   */
    s->t += s->dt;
    s->step += 1;

    MPI_Barrier (MPI_COMM_WORLD);

    return;
}
