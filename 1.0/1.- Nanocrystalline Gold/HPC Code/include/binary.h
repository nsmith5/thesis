#include <fftw3.h>
#include <signal.h>
#include <gsl/gsl_rng.h>
#include <hdf5.h>

/* State structure and Functions */

typedef struct
{
	/* Free Energy Parameters */
	double eta;
	double chi;
	double epsilon0;
	double sigma0;
	double sigma;
	double omega;
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
	ptrdiff_t local_n1;
	ptrdiff_t local_1_start;

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

	/* Random Number Generator */
	gsl_rng *rng;
}
state;

state* create_state (int N, double dx,double dt);
void destroy_state (state* s);
void set_C (state * s);
void make_square (state *s, double h);
void make_const (state *s, double h);
double calc_k (int i, int j, int N, double dx);

/* I/O Functions */

void mpi_print (const char *str);
hid_t io_init_new_file (const char *filename);
hid_t io_init_from_file (const char *filename);
herr_t io_finalize (hid_t file_id);
herr_t save_state (state* s, hid_t file_id);
state* load_state (hid_t file_id, const char *datafile);
state* new_state_from_file (const char *filename);

/* Error handling */

void my_error (const char* error_string);

/* Time stepping Algorithm Functions */

void step (state *s);
void normalize (state *s, fftw_complex *field);

/* Simulation Setup Functions */

extern volatile sig_atomic_t time_to_leave;
void init (int argc, char **argv);
void finalize (void);
