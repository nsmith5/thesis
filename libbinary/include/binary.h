#include <fftw3.h>
#include <signal.h>
#include <gsl/gsl_rng.h>
#include <hdf5.h>

/* State structure and Functions */

/* Field with real and Fourier representations */
typedef struct
{
	double* real;
	fftw_complex* fourier;

}
field;

field* create_field (ptrdiff_t);
void destroy_field (field*);

void fft (fftw_plan, field*);
void ifft (fftw_plan, field*);

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
	double dx;					// Grid spacing
	double dt;					// Time step size
	double t;					// "time"
	int step;					// Time step number
	ptrdiff_t N;				// Size of domain (N x N)
	ptrdiff_t local_n0;			// Local number of rows in real space
	ptrdiff_t local_0_start;	// Local start row in real space
	ptrdiff_t local_n1;			// Local number of rows in fourier space
	ptrdiff_t local_1_start;	// Local start row in fourier space

	/* Fields */
	field *c;		// Concentration
	field *n;		// Reduced density (rho - rho0) / rho
	
	field *cnl;	 	// Nonlinear part of concentration EOM
	field *nnl;  	// Nonlinear part of density EOM
	field *Cn;	 	// C * n term

	fftw_complex *fxin;		// Fourier space noise for density
	fftw_complex *fxic;		// Fourier space noise for concentration

	/* Propagator fields */
	double *Pn;
	double *Qn;
	double *Ln;
	double *Pc;
	double *Qc;
	double *Lc;

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

state* create_state (int, double, double);
void destroy_state (state*s);

double calc_k(int, int, int, double);
void set_C (state*);
void set_k2 (state*);
void set_propagators (state*);

void set_uniform (double*, double, ptrdiff_t, ptrdiff_t);

/* I/O Functions */

void mpi_print (const char*);
hid_t io_init_new_file (const char*);
hid_t io_init_from_file (const char*);
herr_t io_finalize (hid_t);
herr_t save_state (state*, hid_t);
state* load_state (hid_t, const char *);
state* new_state_from_file (const char*);

/* Error handling */

void my_error (const char*);

/* Time stepping Algorithm Functions */

void step (state *s);

/* Simulation Setup Functions */

extern double start_time;
void init (int argc, char **argv);
void finalize (void);
