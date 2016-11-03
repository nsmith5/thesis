#include <fftw3.h>

/* Definition of the simulation state */
typedef struct
{
	/* Free Energy Parameters */
	double eta;
	double chi;
	double epsilon0;
	double sigma0;
	double sigma;
	double omega;
	double kbT;
	double Wc;

	/* Dynamic Parameters */
	double Mn;
	double Mc;

	/* Correlation Function Parameters */
	double k0;
	double alpha;
	double beta;
	double rho;
	double alphac;

	/* Numerical Parameters */
	double dx;
	double dt;
	double t;
	int step;
	ptrdiff_t N;
	ptrdiff_t local_n0;
	ptrdiff_t local_0_start;

	/* Fields */
	double *c;
	double *n;

	/* Temporary fields */
	fftw_complex *fnnl;
	fftw_complex *fcnl;
	fftw_complex *fc;
	fftw_complex *fn;
	fftw_complex *fCn;
	fftw_complex *fxin;
	fftw_complex *fxic;

	double *Cn;
	double *nnl;
	double *cnl;

	/* Operators */
	double *k2;
	double *C;

	/* Fourier transform plans */
	fftw_plan fft_plan;
	fftw_plan ifft_plan;
}
state;

state* create_state (int N, double dx,double dt);

void destroy_state (state* s);

void set_C (state * s);

void make_square (state *s, double h);

void make_const (state *s, double h);

double calc_k (int i, int j, int N, double dx);
