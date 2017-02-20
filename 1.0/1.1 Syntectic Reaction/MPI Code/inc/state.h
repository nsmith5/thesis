
typedef struct {
    // Free Energy Parameters
    double eta;
    double chi;
    double Wc;
    double eps0;
    double sig0;
    double c0;
    double sigma;
    double omega;
    double kbT;
    
    // Dynamic Parameters
    double Mn;
    double Mc;
    
    // Correlation Function Parameters
    double k0;
    double alpha;
    double rho;
    double alpha_c;
    
    // Numerical Parameters
    double dx;
    double dt;
    int N;
    
    // Temporary Fields
    complex* fnnl;
    complex* fcnl;
    complex* fc;
    complex* fn;
    complex* kCn;
    complex* zeta_n;
    complex* zeta_

}
state;
