#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <getopt.h>
#include <random>
#include <thread>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>
#include <string>
#include <chrono>

#include "matrix.h"

#include <boost/math/distributions/exponential.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/random.hpp>

using namespace std;
using boost::math::exponential;
using ldouble = long double;

const int MAXBINS = 1000;
const int MAXERROR = 1e-20;
const int MAXWINDOWSIZE = 1e6;


/* flags */
static int makeplot;
static int nthreads; // number of threads
static int MAXTHREADS; // number of max threads; has be less than the number of bins & CPU currency
static int lvals_per_thread; // number of lvals each thread handles for first (nthread-1) threads
static int lvals_remainder; // remaining lvals are run at final thread call
static int fast_flag;
static int ind_flag; // non-independence assumption
static int mult_flag; // calculate with multiple reference
static int fixed_window_size; // instead of calculating unconditional P to length, pick a fixed length
static int bin_flag;
static int short_flag1; // shortening on branch 1
static int short_flag2; // shortening on branch 2
static int verbose_flag; // if verbose flag is set, then print out the pCorrect values for each L
static int complete_flag; // if complete_flag is on, instead of working with specific params run a complete benchmark
int* thread_args;

string outpath;

/* constants */
// constants for nonind
static ldouble T1 = 2000.0; // this is Tbefore
static ldouble T2 = 10.0; // this is Tafter
static ldouble infty = 100000;
static ldouble t12 = T1+T2;
const ldouble grid_t1 = 100;
const ldouble grid_t2 = 5;
const int ngrid_t1 = (t12-T2)/grid_t1;
const int ngrid_t2 = infty/grid_t2;

// for running with complete_flag
ldouble min_alpha = 0.01;
ldouble alpha_delta = 0.05;
ldouble max_alpha = 0.11;
ldouble min_ne1 = 1000.0;
ldouble max_ne1 = 16000.0;
ldouble ne1_delta = 5000.0;
ldouble min_r = pow(10.0, -10);
ldouble max_r = pow(10.0, -8);
ldouble r_delta = 0.25;
ldouble min_u = pow(10.0, -10);
ldouble max_u = pow(10.0, -8);
ldouble u_delta = 0.25;
ldouble min_tbefore = 1000.0;
ldouble max_tbefore = 21000.0;
ldouble min_tafter = 500.0;
ldouble max_tafter = 10500.0;
ldouble min_tshort = 0.0;
ldouble max_tshort = 5000.0;
ldouble tshort_delta = 100.0;
ldouble intervals = 20.0;
int alpha_ind = 1;
int ne1_ind = 1;
int r_ind = 1;
int u_ind = 1;

//branch shortening; this is the time that the test sample is longer than ref sample
ldouble T_short[2] = {.0, .0};

// constants for ind
static ldouble Tbefore = 12000.0;
static ldouble Tafter = 1000.0;
const ldouble Ne = 10000.0;
ldouble u = 1.25*pow(10.0, -8);
ldouble r = pow(10.0, -8);
ldouble alpha = 0.5;

int t1_max = 1;
int t2_max = 1;
int n_max = 3;
int s_max = 50+1;

/* global variables */
static int numRef=1;
int numBin;
int ct1; // current values
int ct2;

ldouble cTbefore; //current values
ldouble cTafter[2] = {Tafter, Tafter}; // for poisson calculations
ldouble cTafter_real[2] = {Tafter, Tafter}; 
ldouble cTafter_test = Tafter; // for length calculations

/* variables for extrapolation */
int* s_ext_0;

/* for very large tau values, ignore first s til machine precision */
int* s_start_match;
int* s_start_miss;
int* s_ext_match;
int* s_ext_miss;
ldouble* PMatched_rates; // rates of pmatch
ldouble* PMissed_rates;


ldouble ratio_threshold = 0.05;
ldouble ratio_start_threshold = 0.0001;
const ldouble ratio_max_threshold = 0.5;
const int tau_threshold = 5000; // for very large tau, we ignore first s values where gamma ~ 0

int offset_l; // offset in l-indexing to save vals for Population 2
int max_length_pop1; // if maximum of pop1 segment length is smaller than an L of pop2, we just assign it to pop2

int max_s_ext_max; // dynamic alloc for s until max_s_ext_max, which is max(s_ext_max)
const int min_s_ext = 50;
const int max_s = 1000000; // might have to change so it's based on L values
ldouble* eta_0;
ldouble* tau_0;
ldouble* C_0;
ldouble* rho_0;
ldouble* psum_12;
ldouble* psum_11;

ldouble* eta_1;
ldouble* tau_1;
ldouble* t_1;
ldouble* C_1;


ldouble Ne0 = Ne;
ldouble Ne1 = Ne; // assume eff. pop size all same; can modify easily
ldouble Ne2 = Ne;
ldouble* Ne_sizes; // length of l
ldouble* Ne_vals; // length of 2

ldouble* L1;
ldouble* L2;
ldouble* p; //store CDF P(L_i) ~ exp for L_i
ldouble* L;
ldouble* S; // v4.1: now only until s = s_ext+1

/* store the threads here */
thread* mythreads;

/* arrays for comp; multi-thread ver. */
ldouble*** Gamma1;
ldouble*** Gamma2;
ldouble** C1;
ldouble** C2;

ldouble*** Gamma;
ldouble** C;

ldouble** PMatched;

ldouble** PMissMatched;

ldouble*** Result;

ldouble** SumRes;

/* variables for non-independence */
ldouble** poiss1;
ldouble** poiss2;
Vector pCor;

/* Helper functions for array manipulation & I/O */
bool check_file(const string& name){
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}


void printList(ldouble* ptr, int len){
    for (int i=0; i<len; i++)
        cout << fixed << setprecision(5) << ptr[i] << " ";
    cout << endl;
}

void printMatrix(ldouble** ptr, int nrows, int ncols, bool ind, string header){
    for (int i=0; i<nrows; i++){
        if (ind)
            cout << (i+1)*(int)Tafter << "\t";
        for (int j=0; j<ncols; j++)
            cout << header << fixed << setprecision(5) << ptr[i][j] << " ";
        cout << endl;
    }
    cout << endl;
}

void outputMatrix(ldouble** ptr, int nrows, int ncols, string outpath, string header, bool ishead){
    ofstream outfile;
    outfile.open(outpath, ios::out | ios::app);
    if (ishead)
        outfile << "alpha,Ne1,r,u,Tbefore,Tafter,Tshort,pCorrect"<< endl;
    //outfile << "T1: " << cTbefore << " T2, pop1: " << cTafter_real[0] << " T2, pop2: " << cTafter_real[1] << " T2, test: " << cTafter_test << endl;
    for (int i=0; i<nrows; i++){
        for (int j=0; j<ncols; j++)
            outfile << header << to_string(int(Tafter*(i+1))) << "," << 0 << "," << fixed << setprecision(6) << ptr[i][j] << " ";
        outfile << endl;
    }
}

void outputList(ldouble* ptr, int len, string outpath){
    ofstream outfile;
    //outfile.open(outpath, ios::out | ios::app);
    outfile.open(outpath, ios::out);
    for (int i=0; i<len; i++)
        outfile << fixed << setprecision(10) << ptr[i] << " ";
    outfile << endl;
}
void clearList(ldouble* list, int len){
    for (int i=0; i<len; i++)
        list[i] = .0;
}

int argparse(int argc, char* argv[]){
    // argument parsing
    int c;
    if (argc > 1){
        cout << "Options: ";
        for (int i=1; i<argc; i++)
            cout << argv[i] << " ";
        cout << endl;
    }
    while (1){
        static struct option long_options[] = {
            {"plot", no_argument, &makeplot, 1},
            {"branch-short-1", required_argument, &short_flag1, 's'},
            {"branch-short-2", required_argument, &short_flag2, 'S'},
            {"out", required_argument, 0, 'o'},
            {"bins", required_argument, 0, 'b'},
            {"window-size", required_argument, 0, 'w'},
            {"threads", required_argument, 0, 't'},
            {"mutations", required_argument, 0, 'm'},
            {"num-ref", required_argument, 0, 'n'},
            {"tafter", required_argument, 0, 'A'},
            {"tbefore", required_argument, 0, 'B'},
            {"alpha", required_argument, 0, 'a'},
            {"pop-eff", required_argument, 0, 'e'},
            {"pop1", required_argument, 0, 'p'},
            {"pop2", required_argument, 0, 'P'},
            {"mut-rate", required_argument, 0, 'u'},
            {"rec-rate", required_argument, 0, 'r'},
            {"cond", no_argument, &fast_flag, 1},
            {"non-ind", no_argument, &ind_flag, 1},
            {"mult-ref", no_argument, &mult_flag, 1},
            {"verbose", no_argument, &verbose_flag, 1},
            {"complete", required_argument, 0, 'c'},
            {0, 0, 0, 0}
        };
        int option_idx = 0;
        c = getopt_long(argc, argv, "s:S:o:b:w:t:n:m:a:A:B:e:p:P:v:r:u:c:", long_options, &option_idx);
        if (c==-1)
            break;
        switch(c){
            case 0:
                if (long_options[option_idx].flag != 0)
                    break;
            case 'c':
                intervals = atof(optarg);
                cout << "****running on different parameters with " << intervals << " intervals****" << endl;
                complete_flag = 1;
                t1_max= (int)intervals;
                t2_max = (int)intervals;
                break;
            case 's':
                T_short[0] =atof(optarg);
                cout << "branch shortening on pop1: " << T_short[0] << endl;
                if ((T_short[0] < 0) || (T_short[0] >= Tafter)){
                    cerr << "branch shortening time must be positive and smaller than Tafter!" << endl;
                    return EXIT_FAILURE;
                }
                break;
            case 'S':
                T_short[1] =atof(optarg);
                cout << "branch shortening on pop2: " << T_short[1] << endl;
                if ((T_short[1] < 0) || (T_short[1] >= Tafter)){
                    cerr << "branch shortening time must be positive and smaller than Tafter!" << endl;
                    return EXIT_FAILURE;
                }
                break;
            case 'o':
                // outfile to print out the result
                outpath = optarg;
                if ( outpath == "" ){
                    cerr << "output path cannot be empty!" << endl;
                    return EXIT_FAILURE;
                }
                if (check_file(outpath)){
                    cout << "output file already exists!" << endl;
                    return EXIT_FAILURE;
                }
                break;
            case 'b':
                numBin = atoi(optarg)+1;
                if ((numBin < 2) || (numBin > MAXBINS)){
                    cerr << "invalid number of bins!" << endl;
                    return EXIT_FAILURE;
                }
                else if (fixed_window_size){
                    cerr << "cannot have fixed length (conditional) and bin (unconditional) at the same time!"
                    << endl;
                    return EXIT_FAILURE;
                }
                bin_flag = 1;
                break;
            case 'w':
                fixed_window_size = atoi(optarg);
                if ((fixed_window_size < 2) || (numBin > MAXWINDOWSIZE)){
                    cerr << "invalid window size!" << endl;
                    return EXIT_FAILURE;
                }
                else if (numBin){
                    cerr << "cannot have fixed length (conditional) and bin (unconditional) at the same time!"
                    << endl;
                    return EXIT_FAILURE;
                }
                else if (nthreads > 1){
                    cerr << "cannot have multithreading for fixed length!" << endl;
                    return EXIT_FAILURE;
                }
                break;
            case 'm':
                s_max = atoi(optarg)+1;
                cout << "max # of mutations: " << s_max-1 << endl;
                if ((s_max < min_s_ext) || (s_max > max_s)){
                    cerr << "invalid number of mutations!" << endl;
                    return EXIT_FAILURE;
                }
                break;
            case 'n':
                numRef = atoi(optarg);
                cout << "unit number of references: " << numRef << endl;
                if ((numRef < 1) || (numRef > 1000)){
                    cerr << "invalid number of references! Must be an integer between [1, 100]" << endl;
                    return EXIT_FAILURE;
                }
                break;
            case 'A':
                Tafter = atoi(optarg);
                cout << "max Tafter: " << Tafter*t1_max << endl;
                if (Tafter < 0){
                    cerr << "invalid Tafter value!" << endl;
                    return EXIT_FAILURE;
                }
                break;
            case 'B':
                Tbefore = atoi(optarg);
                cout << "max Tbefore: " << Tbefore*t2_max << endl;
                if (Tbefore < 0){
                    cerr << "invalid Tbefore value!" << endl;
                    return EXIT_FAILURE;
                }
                break;
            case 't':
                nthreads = atoi(optarg);
                MAXTHREADS = thread::hardware_concurrency();
                // if arg is 0, then just set optimal # of threads
                if (!nthreads){
                    cout << "automatically setting # of threads to: " << MAXTHREADS <<  endl;
                    nthreads = MAXTHREADS;
                }
                if ((nthreads < 1) || (nthreads > MAXTHREADS)){
                    cerr << "invalid number of threads! (must be greater than 0 and less than " << MAXTHREADS << endl;
                    return EXIT_FAILURE;
                }
                if (fixed_window_size){
                    cerr << "cannot have multithreading for fixed length!" << endl;
                    return EXIT_FAILURE;
                }
                break;
            case 'a':
                alpha = atof(optarg);
                if ((alpha >= 1.0) || (alpha < 0.0)){
                    cerr << "invalid alpha value! (must be between 0 and 1.0)" << endl;
                    return EXIT_FAILURE;
                }
                cout << "alpha: " << alpha << endl;
                break;
            case 'e':
                Ne0 = atof(optarg);
                if ((Ne0 < 100) || (Ne0 > 50000.0)){
                    cerr << "invalid Ne0 (effective population size 0) value! (must be between 100 and 50,000)" << endl;
                    return EXIT_FAILURE;
                }
                cout << "Ne0: " << Ne0 << endl;
                break;
            case 'p':
                Ne1 = atof(optarg);
                if ((Ne1 < 100) || (Ne1 > 50000.0)){
                    cerr << "invalid Ne1 (effective population size 1) value! (must be between 100 and 50,000)" << endl;
                    return EXIT_FAILURE;
                }
                cout << "Ne1: " << Ne1 << endl;
                break;
            case 'P':
                Ne2 = atof(optarg);
                if ((Ne2 < 100) || (Ne2 > 50000.0)){
                    cerr << "invalid Ne2 (effective population size 2) value! (must be between 100 and 50,000)" << endl;
                    return EXIT_FAILURE;
                }
                cout << "Ne2: " << Ne2 << endl;
                break;
            case 'u':
                u = atof(optarg);
                if ((u < 0) || (u > 1e-5)){
                    cerr << "invalid mutation rate value! (must be between 0 and 1e-5)" << endl;
                    return EXIT_FAILURE;
                }
                min_u = u;
                cout << "u: " << u << endl;
                break;
            case 'r':
                r = atof(optarg);
                if ((r < 0) || (r > 1e-5)){
                    cerr << "invalid recombination rate value! (must be between 0 and 1e-5)" << endl;
                    return EXIT_FAILURE;
                }
                cout << "r: " << r << endl;
                break;
            default:
                abort();
        }
    }
    if (fast_flag)
        cout << "calculating P conditional to the smaller population" << endl;
    if (ind_flag)
        cout << "calculating without the independence of lineage assumption" << endl;
    if (mult_flag)
        cout << "calculating for multiple reference panels" << endl;
    if (mult_flag && ind_flag){
        cerr << "Cannot have independence & multiple references flags at the same time" << endl;
        return EXIT_FAILURE;
    }

    if ((nthreads >= numBin) && (nthreads)){
        cerr << "Too many threads; nthread must be smaller than numBin" << endl;
        return EXIT_FAILURE;
    }
    if (!nthreads){
        cout << "Running on single thread" << endl;
        nthreads = 1;
    }
    if (numBin){
        if (nthreads != 1){
            lvals_per_thread = ceil(float(numBin-1)/nthreads);
            lvals_remainder = lvals_per_thread - (lvals_per_thread*nthreads - (numBin-1));
        // if last thread has less than half of other threads
            if (lvals_remainder < ceil(float(lvals_per_thread)/2.0)){
                lvals_per_thread -= 1;
                lvals_remainder = (numBin-1) - (lvals_per_thread*(nthreads-1));
            }
        }
        else
            lvals_remainder = numBin-1;
        cout << "number of bins:  " << numBin - 1 << endl << "number of threads: " << nthreads << endl;
        cout << "bins per thread (first " << nthreads-1  << " threads): " << lvals_per_thread << endl 
            << "bins remaining (last thread): " << lvals_remainder << endl;

    }
    else{
        lvals_per_thread = 0;
        lvals_remainder = 1;
    }
    // assign thread args
    thread_args = new int[nthreads];
    for (int i=0; i < nthreads; i++)
        thread_args[i] = i;
    cout << "Parsed all options." << endl;
    return 0;
}

ldouble getSum(ldouble* list, int len){
    ldouble sum = .0;
    for(int i=0; i<len; i++)
        sum += list[i];
    return sum;
}

ldouble getMin(ldouble* l1, ldouble* l2, int len){
    ldouble minval = 99999999.9;
    for (int i=0; i<len; i++)
        minval = min(min(l1[i], l2[i]), minval);
    return minval;
}

ldouble getMax(ldouble* l1, ldouble* l2, int len){
    ldouble maxval = -99999999.9;
    for (int i=0; i<len; i++)
        maxval = max(max(l1[i], l2[i]), maxval);
    return maxval;
}

// expinv function in MATLAB; fills up result array with icdf for pvals
void expinv(ldouble* pvals, int npvals, ldouble mu, ldouble* result){
    exponential exp(1/mu);
    for(int i=0; i<npvals; i++)
        result[i] = quantile(exp, pvals[i]);
}

void allocateMemExt(int length){
    s_start_match = new int[length];
    s_start_miss = new int[length];
    s_ext_match = new int[length];
    s_ext_miss = new int[length];
    s_ext_0 = new int[length];
    eta_0 = new ldouble[length];
    tau_0 = new ldouble[length];
    C_0 = new ldouble[length];
    rho_0 = new ldouble[length];
    psum_12 = new ldouble[length];
    psum_11 = new ldouble[length];
    eta_1 = new ldouble[length];
    tau_1 = new ldouble[length];
    C_1 = new ldouble[length];
    t_1 = new ldouble[length];
    Ne_sizes = new ldouble[length];
    Ne_vals = new ldouble[2];
    Ne_vals[0] = Ne1;
    Ne_vals[1] = Ne2;
    // init
    for (int l=0; l < length; l++){
        s_start_match[l] = 0;
        s_start_miss[l] = 0;
        s_ext_match[l] = 0;
        s_ext_miss[l] = 0;
        s_ext_0[l] = .0;
        eta_0[l] = .0;
        tau_0[l] = .0;
        C_0[l] = .0;
        rho_0[l] = .0;
        psum_12[l] = .0;
        psum_11[l] = .0;
        eta_1[l] = .0;
        tau_1[l] = .0;
        C_1[l] = .0;
        t_1[l] = .0;
        Ne_sizes[l] = (l < offset_l) ? Ne1 : Ne2;
    }
}

void cleanUpMemExt(){
    delete[] s_start_match;
    delete[] s_start_miss;
    delete[] s_ext_match;
    delete[] s_ext_miss;
    delete[] s_ext_0;
    delete[] eta_0;
    delete[] tau_0;
    delete[] C_0;
    delete[] rho_0;
    delete[] psum_12;
    delete[] psum_11;
    delete[] eta_1;
    delete[] tau_1;
    delete[] C_1;
    delete[] t_1;
    delete[] Ne_sizes;
    delete[] Ne_vals;
}

/* calculate length distribution for given params & allocate length memories */
void setLengthDistribution(ldouble cTa){
    ldouble mean1 = 1/(r*(1-alpha)*cTa);
    ldouble mean2 = 1/(r*alpha*cTa);
    if (!numBin){
        if (!fixed_window_size){
            // simply calculate P(corr|L) where L is s.t. P(L) = 0.5
            p = new ldouble[2];
            p[0] = 0.5;
            p[1] = 0.975;
            L1 = new ldouble[2];
            L2 = new ldouble[2];
            L = new ldouble[4];
            expinv(p, 2, mean1, L1);
            expinv(p, 2, mean2, L2);
            offset_l = 1;
            for (int i=0; i < 2; i++){
                L[i] = L1[i];
                L[i+2] = L2[i];
            }
            numBin = 2;
        }
        // fixed window (conditional to fixed length)
        else{
            // for simplicity of code, we put L[1] & L[3] dummy vals
            // so it's the same as numBin==2 above
            L = new ldouble [4];
            for (int i=0; i<4; i++)
                L[i] = (ldouble)fixed_window_size;
            offset_l = 1;
            numBin = 2;
        }
    }
    else{
        offset_l = numBin-1;
        p = new ldouble[numBin-1];
        for(int i=1; i<numBin; i++)
            p[i-1] = (ldouble)i/numBin;
        L1 = new ldouble[numBin-1];
        L2 = new ldouble[numBin-1];
        L = new ldouble[numBin-1+offset_l];
        expinv(p, numBin-1, mean1, L1);
        expinv(p, numBin-1, mean2, L2);
        for (int i=0; i<numBin-1; i++){
            L[i] = L1[i];
            L[i+offset_l] = L2[i];
        }
    }
    // v4: allocate memory for extrapolation variables
    if (!ind_flag)
        allocateMemExt(numBin-1+offset_l);
}

/* for different Tafter, length distribution changes */
void updateLengthDistribution(double cTa, bool verbose){
    ldouble mean1 = 1/(r*(1-alpha)*cTa);
    ldouble mean2 = 1/(r*alpha*cTa);
    if (!bin_flag){
        if (!fixed_window_size){
            expinv(p, 2, mean1, L1);
            expinv(p, 2, mean2, L2);
            for (int i=0; i < 2; i++){
                L[i] = L1[i];
                L[i+2] = L2[i];
            }
        }
    }
    else{
        expinv(p, numBin-1, mean1, L1);
        expinv(p, numBin-1, mean2, L2);
        for (int i=0; i<numBin-1; i++){
            L[i] = L1[i];
            L[i+offset_l] = L2[i];
            if (verbose)
                cout << L[i] << " " << L[i+offset_l] << endl;
        }
    }
}

// gammainc function in MATLAB; fills up result array with upper incomplete gamma
void ugammainc(ldouble rate, int slen, ldouble* result, int offset){
    for(int i=0; i<slen; i++)
        result[i] = boost::math::gamma_q(i+1.0+offset, rate);
}

// gammainc function in MATLAB; fills up result array with lower incomplete gamma
void lgammainc(ldouble rate, int slen, ldouble* result, int offset){
    for(int i=0; i<slen; i++)
        result[i] = boost::math::gamma_p(i+1.0+offset, rate);
}

// poisson pdf; fills up result array with poisson pdf values for different s values
void poisson_pdf(ldouble rate, int slen, ldouble* result, int offset){
    boost::math::poisson_distribution<> poisson(rate);
    for(int i=0; i<slen; i++)
        result[i] = boost::math::pdf(poisson, i+offset);
}

// for large rate, ignore first s values
int findStartPoint(ldouble rate){
    int starts = (int) rate;
    if (starts > 500){
        ratio_threshold = min(3.168631975852/rate-795.2743812/pow(rate, 2)+0.0001000489, ratio_max_threshold);
        for (int s=starts; s > -1; s--){
            ldouble denom = abs(log(starts) - log(s) - 1/(2*s)+1/(8*pow(s,2)));
            ldouble norm_ratio = 1/denom/(2*starts);
            if (norm_ratio <  ratio_start_threshold){
                if (boost::math::gamma_q(s+1.0, rate) < 1e-24){
                    return s;
                    break;
                }
            }
        }

    }
    return 0;
}

int findExtPoint(ldouble rate, int offset){
    // v6: select s_ext analytically: we numerically calculate when is the good threshold point for some sampled values of tau,
    // fit a curve and use that as a threshold (instead of some static value across all tau, which works poorly)
    int starts = (int)rate;
    if (starts > 1){
        ratio_threshold = min(0.3966873131/rate-0.69712127028/pow(rate, 2)+0.527098880906/pow(rate, 3)+0.000706932491, ratio_max_threshold);
        for (int s=starts; s < s_max+offset; s++){
            ldouble denom = abs(log(starts) - log(s) - 1/(2*s)+1/(8*pow(s,2)));
            ldouble norm_ratio = 1/denom/(2*starts);
            if (norm_ratio < ratio_threshold) {
                if (1.0 - boost::math::gamma_q(s+1.0, rate) < 1e-24){
                    if (s < s_max){
                        return max(s, min_s_ext);
                    }
                }
            }
        }
    }
    else
        return min_s_ext;
    return s_max+offset;
}

// set global variables for a given Tafter and Tbefore values
void setGlobalVars(){
    if (!mult_flag){
        int lmax = (fast_flag) ? numBin-1 : numBin-1+offset_l;
        Ne_vals[0] = Ne1;
        Ne_vals[1] = Ne2;
        for(int l=0; l <lmax; l++){
            Ne_sizes[l] = (l < offset_l) ? Ne1 : Ne2;
            int pop = l/offset_l;
            eta_0[l] = 4*Ne0*L[l]*u+1;
            tau_0[l] = (cTbefore+cTafter[pop])*(1/(2*Ne0) + 2*L[l]*u);
            C_0[l] = exp((cTafter[pop]+cTbefore)/(2*Ne0))/eta_0[l];
            rho_0[l] = (eta_0[l] - 1)/eta_0[l];
            
            eta_1[l] = 4*Ne_vals[l/offset_l]*L[l]*u+1;
            tau_1[l] = (cTbefore+cTafter[pop])*(1/(2*Ne_vals[l/offset_l]) + 2*L[l]*u);
            t_1[l] = cTafter[pop]*(1/(2*Ne_vals[l/offset_l]) + 2*L[l]*u);
            C_1[l] = exp((cTbefore+cTafter[pop])/(2*Ne0) - cTbefore/(2*Ne_vals[l/offset_l]))/eta_0[l];
            s_ext_0[l] = min(findExtPoint(max(max(tau_0[l],tau_1[l]), t_1[l]),0 ), s_max-1);
        }
    }
    max_length_pop1 = getMax(L, L, offset_l);
    if (verbose_flag)
        cout << "Max length_pop1: " << max_length_pop1 << endl;
}

/* calculate max_s_ext_max, in order to allocate array memory */
void setMaxS(){
    if (mult_flag || complete_flag){
        max_s_ext_max = s_max+1000;
    }
    else{
        int maxval = 0;
        int lmax = (fast_flag) ? numBin-1 : numBin-1+offset_l;
        int t2 = t2_max;
        for (int t1 = 1; t1 < t1_max+1; t1++){
            int Tb = Tbefore*t2;
            int Ta1 = Tafter*t1 - T_short[0]/2;
            int Ta2 = Tafter*t1 - T_short[1]/2;
            int Ta = Tafter*t1;
            updateLengthDistribution(Ta, false);
            for (int l=0; l < lmax; l++){
                int rate1 = l < offset_l ? (Tb+Ta1)*(1/(2*Ne0) + 2*L[l]*u) : (Tb+Ta2)*(1/(2*Ne0) + 2*L[l]*u);
                int rate2 = l < offset_l ? (Tb+Ta1)*(1/(2*Ne1) + 2*L[l]*u) : (Tb+Ta2)*(1/(2*Ne2) + 2*L[l]*u);
                int rate3 = l < offset_l ? (Ta1)*(1/(2*Ne1) + 2*L[l]*u) : (Ta2)*(1/(2*Ne2) + 2*L[l]*u);
                maxval = max(findExtPoint(max(max(rate1, rate2), rate3), 0), maxval);
            }
        }
        max_s_ext_max = maxval+1000;
    }
    cout << "max_s_ext: " << max_s_ext_max << endl;
}

// anc == 0 --> pop1, anc == 1 --> pop2
ldouble getExpectedRefLineage(int numref, bool anc, bool ref){
    if (anc == ref)
        return 4*numref*Ne_vals[ref]/(numref*cTafter_real[ref] + 4*Ne_vals[ref]);
    else
        return 4*numref*Ne_vals[ref]/(numref*(cTafter_real[ref]+cTbefore) + 4*Ne_vals[ref]);
}

ldouble getExpectedCoalesceT(int numref, bool anc, bool ref){
    if (anc == ref)
        return cTafter_real[ref]+ 4*Ne_vals[ref]/(getExpectedRefLineage(numref, anc, ref)+1);
    else
        return cTafter_real[ref]+cTbefore+4*Ne0/(getExpectedRefLineage(numref, anc, ref)+1);
}

void getProbMatchedRate_mult(int numref, int l, int tid, bool anc){
    ldouble L_val = L[l];
    ldouble rate = 2*u*L_val*getExpectedCoalesceT(numref, anc, anc);
    s_start_match[l] = findStartPoint(rate);
    s_ext_match[l] = findExtPoint(rate, s_start_match[l]);
    PMatched_rates[tid] = rate;
    int range = s_ext_match[l] - s_start_match[l];
    if (range > s_max)
        cout << "need more mutations to sum over: " << range << endl;
}

void getProbLargerMissMatchedRefPanelRate_mult(int numref, int l, int tid, bool anc){
    ldouble L_val = L[l];
    ldouble rate = 2*u*L_val*getExpectedCoalesceT(numref, anc, !anc);
    s_start_miss[l] = findStartPoint(rate);
    s_ext_miss[l] = findExtPoint(rate, s_start_match[l]);
    PMissed_rates[tid] = rate;
    int range = s_ext_miss[l] - s_start_miss[l]; 
    if (range > s_max)
        cout << "need more mutations to sum over: " << range << endl;
}

void getProbMatched_mult(int l, int tid){
    poisson_pdf(PMatched_rates[tid], s_max, PMatched[tid], min(s_start_match[l], s_start_miss[l]));
}

void getProbLargerMissMatchedRefPanel_mult(int l, int tid){
    lgammainc(PMissed_rates[tid], s_max, PMissMatched[tid], min(s_start_match[l], s_start_miss[l]));
    int range = max(s_ext_miss[l], s_ext_match[l]) - min(s_start_match[l], s_start_miss[l]); 
    if (range > s_max)
        cout << "need more mutations to sum over: " << range << endl;
}


/*** single reference panel calculations ***/
/* calculate PMiss with extrapolation*/
ldouble extrapolatePMiss(int lenS, int l){
    ldouble sumPMissMatched = psum_12[l] + C_0[l]*((pow(rho_0[l], lenS) - pow(rho_0[l], s_ext_0[l]))/(rho_0[l]-1) + 0.5*pow(rho_0[l], lenS));
    return (1 - sumPMissMatched);
}

ldouble extrapolatePMatch(int lenS, int l){
    ldouble PMatchedval = C_1[l]*pow(rho_0[l], lenS);
    return PMatchedval;
}

void getProbMatched(int lenS, int l, int tid, int pop){
    ldouble L_val = L[l];
    C1[tid][lenS] = exp(cTafter[pop]/(2*Ne_sizes[l]))*pow(((4*Ne_sizes[l]*L_val*u)/(1+4*L_val*u*Ne_sizes[l])),S[lenS])/(1+4*L_val*u*Ne_sizes[l]);
    C2[tid][lenS] = exp((cTbefore+cTafter[pop])/(2*Ne0)-cTbefore/(2*Ne_sizes[l]))*pow((4*L_val*u*Ne0)/(1+4*L_val*u*Ne0),S[lenS])/(1+4*L_val*u*Ne0);
    PMatched[tid][lenS] = C1[tid][lenS]*(Gamma1[tid][pop][lenS]-Gamma2[tid][pop][lenS]) + C2[tid][lenS]*Gamma[tid][pop][lenS];
}

ldouble getProbMatchedRefPanel(int s, int l, int tid, int pop){
    getProbMatched(s, l, tid, pop);
    if (s == s_ext_0[l])
        psum_11[l] = getSum(PMatched[tid], s);
    return PMatched[tid][s];
}

void getProbMissMatched(int lenS, int l, int tid, int pop){
    ldouble L_val = L[l];
    C[tid][lenS] = exp((cTbefore+cTafter[pop])/(2*Ne0))*pow(((4*Ne0*L_val*u)/(1+4*L_val*u*Ne0)), S[lenS])/(1+4*L_val*u*Ne0);
    PMissMatched[tid][lenS] = C[tid][lenS]*Gamma[tid][pop][lenS];
}

ldouble getProbLargerMissMatchedRefPanel(int s, int l, int tid, int pop){
    getProbMissMatched(s, l, tid, pop);
    ldouble PMissAllLarger = 1-getSum(PMissMatched[tid], s)-0.5*PMissMatched[tid][s];
    if (s == s_ext_0[l])
        psum_12[l] = getSum(PMissMatched[tid], s);
    return PMissAllLarger;
}

// allocate Gamma values at once
void getGammaVals(int l, int tid){
    ldouble L_val1 = L[l];
    ldouble L_val2 = L[l+offset_l];

    ldouble Rate10 = (cTbefore+cTafter[0])*(1/(2*Ne0) + 2*L_val1*u);
    ldouble Rate11 = cTafter[0]*(1/(2*Ne1) + 2*L_val1*u);
    ldouble Rate12 = (cTbefore+cTafter[0])*(1/(2*Ne1) + 2*L_val1*u);
    ugammainc(Rate10, min(s_ext_0[l]+100, (int)L[l]), Gamma[tid][0], 0);
    ugammainc(Rate11, min(s_ext_0[l]+100, (int)L[l]), Gamma1[tid][0], 0);
    ugammainc(Rate12, min(s_ext_0[l]+100, (int)L[l]), Gamma2[tid][0], 0);

    if (!fast_flag){
        ldouble Rate20 = (cTbefore+cTafter[1])*(1/(2*Ne0) + 2*L_val2*u);
        ldouble Rate21 = cTafter[1]*(1/(2*Ne2) + 2*L_val2*u);
        ldouble Rate22 = (cTbefore+cTafter[1])*(1/(2*Ne2) + 2*L_val2*u);
        ugammainc(Rate20, min(s_ext_0[l+offset_l]+100, (int)L[l]), Gamma[tid][1], 0);
        ugammainc(Rate21, min(s_ext_0[l+offset_l]+100, (int)L[l]), Gamma1[tid][1], 0);
        ugammainc(Rate22, min(s_ext_0[l+offset_l]+100, (int)L[l]), Gamma2[tid][1], 0);
    }
}

// clear out Result & SumRes matrices for a new t2 value
void clearRes_mt(){
    for (int i=0; i<nthreads; i++){
        for (int j=0; j<t1_max; j++){
            for (int k=0; k<n_max; k++)
                Result[i][j][k] = .0;
        }
    }
    for (int i=0; i<t1_max; i++){
        for (int j=0; j<n_max; j++)
            SumRes[i][j] = .0;
    }
}

// add up all the multiple threads outputs in the Result matrix to SumRes
void summarize(){
    for (int i=0; i<nthreads; i++){
        for (int j=0; j<t1_max; j++){
            for (int k=0; k<n_max; k++)
                SumRes[j][k] += Result[i][j][k];
        }
    }
}

// initialize; allocate memory for all array parameters
void init_mt(){
    // allocate memory only up to max_s_ext_max
    mythreads = new thread[nthreads];
    S = new ldouble[max_s_ext_max];
    for (int i=0; i<max_s_ext_max; i++)
        S[i] = i;
    Gamma1 = new ldouble**[nthreads];
    Gamma2 = new ldouble**[nthreads];
    Gamma = new ldouble**[nthreads];
    C1 = new ldouble*[nthreads];
    C2 = new ldouble*[nthreads];
    C = new ldouble*[nthreads];
    PMatched_rates = new ldouble[nthreads];
    PMissed_rates = new ldouble[nthreads];
    PMissMatched = new ldouble*[nthreads];
    PMatched = new ldouble*[nthreads];
    Result = new ldouble**[nthreads];
    for (int i=0; i<nthreads; i++){
        Gamma1[i] = new ldouble*[2];
        Gamma2[i] = new ldouble*[2];
        Gamma[i] = new ldouble*[2];
        for (int j=0; j<2; j++){
            Gamma1[i][j] = new ldouble[max_s_ext_max];
            Gamma2[i][j] = new ldouble[max_s_ext_max];
            Gamma[i][j] = new ldouble[max_s_ext_max];
            for (int k=0; k<max_s_ext_max; k++){
                Gamma1[i][j][k] = .0;
                Gamma2[i][j][k] = .0;
                Gamma[i][j][k] = .0;
            }
        }
        C1[i] = new ldouble[max_s_ext_max];
        C2[i] = new ldouble[max_s_ext_max];
        C[i] = new ldouble[max_s_ext_max];
        PMissMatched[i] = new ldouble[max_s_ext_max];
        PMatched[i] = new ldouble[max_s_ext_max];
        for (int j=0; j<max_s_ext_max; j++){
            C1[i][j] = .0;
            C2[i][j] = .0;
            C[i][j] = .0;
            PMissMatched[i][j] = .0;
            PMatched[i][j] = .0;
        }

        Result[i] = new ldouble*[t1_max];
        for (int j=0; j<t1_max; j++){
            Result[i][j] = new ldouble[n_max];
            for (int k=0; k<n_max; k++){
                Result[i][j][k] = .0;
            }
        }
    }
    SumRes = new ldouble*[t1_max];
    for (int i=0; i<t1_max; i++){
        SumRes[i] = new ldouble[n_max];
        for (int j=0; j<n_max; j++)
            SumRes[i][j] = .0;
    }
    cout << "initialization complete" << endl;
}

// de-allocate memory at termination
void cleanUp_mt(){
    for (int i=0; i<nthreads; i++){
        for (int j=0; j<2; j++){
            delete[] Gamma1[i][j];
            delete[] Gamma2[i][j];
            delete[] Gamma[i][j];
        }
        delete[] Gamma1[i];
        delete[] Gamma2[i];
        delete[] Gamma[i];
        delete[] C1[i];
        delete[] C2[i];
        delete[] C[i];
        delete[] PMissMatched[i];
        delete[] PMatched[i];
        for (int j=0; j<t1_max; j++){
            delete[] Result[i][j];
        }
        delete[] Result[i];
    }
    delete[] mythreads;
    delete[] thread_args;
    delete[] L;
    delete[] L1;
    delete[] L2;
    delete[] p;
    delete[] S;
    delete[] Gamma1;
    delete[] Gamma2;
    delete[] Gamma;
    delete[] C1;
    delete[] C2;
    delete[] C;
    delete[] PMatched_rates;
    delete[] PMissed_rates;
    delete[] PMissMatched;
    delete[] PMatched;
    delete[] Result;
    for (int i=0; i<t1_max; i++)
        delete[] SumRes[i];
    delete[] SumRes;
    cleanUpMemExt();
    cout << "cleanup complete" << endl;
}

// single reference panel, independence assumption

void* threadfunc_single(int *tidx){
    int tid = *tidx;
    int l0 = tid*lvals_per_thread;
    int l1 = (tid+1)*lvals_per_thread;
    if (tid == nthreads-1)
        l1 = l0 + lvals_remainder;
     
    for (int l=l0; l<l1; l++){
        // gamma only depends on L, so for a given L value, calculate all Gamma values
        getGammaVals(l, tid);
        ldouble pCorrect1 = .0;
        ldouble pCorrect2 = .0;
        for(int s=0; s<min(s_ext_0[l]+1, int(L[l])); s++){
            pCorrect1 += getProbLargerMissMatchedRefPanel(s, l, tid, 0) * getProbMatchedRefPanel(s, l, tid, 0);
        }

        // extrapolation
        for (int s=s_ext_0[l]+1; s < min(s_max, int(L[l])); s++){
            pCorrect1 += extrapolatePMiss(s, l) * extrapolatePMatch(s, l);
        }
        //pCor(l) = pCorrect1;

        // if fast_flag, calculate the conditional on Pop1; otherwise, calculate unconditonal
        if (!fast_flag){
            pCorrect1 *= alpha;
            for (int s=0; s<min(s_ext_0[l+offset_l]+1, int(L[l+offset_l])); s++){
                pCorrect2 += getProbLargerMissMatchedRefPanel(s, l+offset_l, tid, 1)*getProbMatchedRefPanel(s, l+offset_l, tid, 1);
            }
            for (int s=s_ext_0[l+offset_l]+1; s < min(s_max, int(L[l+offset_l])); s++){
                pCorrect2 += extrapolatePMiss(s, l+offset_l) * extrapolatePMatch(s, l+offset_l);
            }

            // combine two probabilities
            pCorrect1 += pCorrect2*(1-alpha);
        }
        Result[tid][ct1-1][0] += pCorrect1/(numBin-1);
    }
}






// multiple references
void* threadfunc_mult(int *tidx){
    int tid = *tidx;
    int l0 = tid*lvals_per_thread;
    int l1 = (tid+1)*lvals_per_thread;
    if (tid == nthreads-1)
        l1 = l0 + lvals_remainder;
    
    for (int l=l0; l<l1; l++){
        for(int n=1; n<n_max+1; n++){
            ldouble pCorrect1 = .0;
            ldouble pCorrect2 = .0;
            getProbMatchedRate_mult(pow(numRef, n), l, tid, 0);
            getProbLargerMissMatchedRefPanelRate_mult(pow(numRef, n), l, tid, 0);
            getProbMatched_mult(l, tid);
            getProbLargerMissMatchedRefPanel_mult(l, tid);
            
            for(int s=0; s<min((ldouble)s_max, L[l]); s++){
                pCorrect1 += PMatched[tid][s]*PMissMatched[tid][s];
            }
      
            pCor(l) = pCorrect1;
            
            if (!fast_flag){
                pCorrect1 *= alpha;

                getProbMatchedRate_mult(pow(numRef, n), l+offset_l, tid, 1);
                getProbLargerMissMatchedRefPanelRate_mult(pow(numRef, n), l+offset_l, tid, 1);
                getProbMatched_mult(l+offset_l, tid);
                getProbLargerMissMatchedRefPanel_mult(l+offset_l, tid);
                
                for(int s=0; s<min((ldouble)s_max, L[l+offset_l]); s++){
                    pCorrect2 +=  PMatched[tid][s]*PMissMatched[tid][s];
                }
                pCorrect1 += pCorrect2*(1-alpha);
            }
            Result[tid][ct1-1][n-1] += pCorrect1/(numBin-1);
        }
    }
}


void init_nonind(){
    poiss1 = new ldouble*[nthreads];
    poiss2 = new ldouble*[nthreads];
    for (int i=0; i<nthreads; i++){
        poiss1[i] = new ldouble[500000];
        poiss2[i] = new ldouble[500000];
    }
    S = new ldouble[500000];
    for (int i=0; i<500000; i++)
        S[i] = i;
    mythreads = new thread[nthreads];
}

void cleanUp_nonind(){
    for (int i=0; i<nthreads; i++){
        delete[] poiss1[i];
        delete[] poiss2[i];
    }
    delete[] poiss1;
    delete[] poiss2;
    delete[] S;
    delete[] p;
    delete[] L1;
    delete[] L2;
    delete[] L;
    delete[] mythreads;
    delete[] thread_args;
}

void* threadfunc_nonind(int *tidx){
    int tid = *tidx;
    int l0 = tid*lvals_per_thread;
    int l1 = (tid+1)*lvals_per_thread;
    if (tid == nthreads-1)
        l1 = l0 + lvals_remainder;
    
    Matrix pMat(ngrid_t2, ngrid_t1);

    for (int l=l0; l<l1; l++){
        int i=0, j=0;
        for (ldouble t2=t12; t2 < t12 + infty; t2+=grid_t2){
            for (ldouble t1=T2; t1 < t12; t1 +=grid_t1){
                ldouble lambda1 = u*L[l]*(t1-T_short[0]);
                ldouble lambda2 = u*L[l]*(2*t2-t1);
                ldouble pcdf = .0;
                int s2max = findExtPoint(lambda2, 0);
                int s1max = findExtPoint(lambda1, 0);
                int offset_s = min(findStartPoint(lambda2), findStartPoint(lambda2));
                poisson_pdf(lambda1, s1max+10-offset_s, poiss1[tid], offset_s);
                poisson_pdf(lambda2, s2max+10-offset_s, poiss2[tid], offset_s);
                for (int s1=0; s1 < min(s1max, (int)L[l])-offset_s; s1+=1){
                    // calculate pois pdf (l2, s2) over s2=s1+1 to \infty
                    // but instead of \infty just use s2max
                    pcdf += poiss1[tid][s1]*(0.5*poiss2[tid][s1]+getSum(poiss2[tid]+s1+1, s2max-s1));
                    //pcdf += poiss1[tid][s1]*(alpha*poiss2[tid][s1]+getSum(poiss2[tid]+s1+1, s2max-s1));
                }
                pcdf *=  exp(-(t1 - T2)/(2*Ne1))*exp(-(t2-T1-T2)/(2*Ne0));
                pMat(i,j) = pcdf;
                j++;
            }
            i++;
            j=0;
        }
        ldouble xscale=10;
        ldouble yscale=10;

        pMat*= grid_t1*grid_t2/(4*Ne1*Ne0);

        ldouble int_res = bicubic_mem_eff(pMat, xscale, yscale);

        pCor(l) = int_res;

        cout << "T1: " << T1 << " T2: " << T2 << " L : " << L[l] << " pcor: " << int_res << ", " << pMat.sum() << endl;

    }
    pMat.output("pmat.txt");
}

/* main routine */
int main(int argc, char *argv[]){
    if (argparse(argc, argv))
        return EXIT_FAILURE;
   
    auto start = std::chrono::high_resolution_clock::now();

    /* testing implementation for non-independence assumption */
    if (ind_flag){
        init_nonind();
        setLengthDistribution(T2);
        
        pCor.resize(numBin-1);
        cout << "pCor shape: (" << pCor.rows() << "," << pCor.cols() << ")" << endl;
        
        for (int tidx=0; tidx < nthreads; tidx++){    
            mythreads[tidx] = thread(&threadfunc_nonind, &thread_args[tidx]);        
        }
        for (int tidx=0; tidx < nthreads; tidx++){
            mythreads[tidx].join();
        }
        if (verbose_flag){
            string Ne_str = to_string((int)Ne1/1000)+'k';
            pCor.output("pCor.nonind."+Ne_str+".txt");
            outputList(L, numBin-1, "Lvals.nonind."+Ne_str+".txt");
        }
        cout << "Tbefore: " << T1 << " Tafter: " << T2 << " pCorrect : " << pCor.sum()/pCor.cols() << endl;

        cleanUp_nonind();
    }
    else if (mult_flag){
        setLengthDistribution(Tafter);
        setMaxS();
        init_mt(); 
        pCor.resize(numBin-1);
        cout << "pCor shape: (" << pCor.rows() << "," << pCor.cols() << ")" << endl;
        
        for(ct2=1; ct2 < t2_max+1; ct2++){
            for(ct1=1; ct1 < t1_max+1; ct1++){
                cTbefore = Tbefore*ct2;
                cTafter_test = Tafter*ct1;
                for (int i=0; i<2; i++){
                    cTafter[i] = Tafter*ct1 - T_short[i]/2;
                    cTafter_real[i] = Tafter*ct1 - T_short[i];
                }
                updateLengthDistribution(cTafter_test, false);
                setGlobalVars();
                for (int tidx=0; tidx < nthreads; tidx++){
                    mythreads[tidx] = thread(&threadfunc_mult, &thread_args[tidx]);
                }
                for (int tidx=0; tidx < nthreads; tidx++){
                    mythreads[tidx].join();
                }
            }
            summarize();
            cout << "Tbefore(T2): " << int(cTbefore) << " Tafter1 (T1): " << int(cTafter_real[0]) << " Tafter2 (T2): " << int(cTafter_real[1]) << endl << "T1\\ref:";
            for (int n=1; n<n_max+1; n++)
                cout << "\t" << pow(numRef, n);
            cout << endl;
            printMatrix(SumRes, t1_max, n_max, true, "");
            if (outpath!="")
                outputMatrix(SumRes, t1_max, n_max, outpath, "", false);
            clearRes_mt();
        }
        if (verbose_flag){
            string Ne_str = to_string((int)Ne1/1000)+'k';
            pCor.output("pCor.ind.mult."+Ne_str+".txt");
            outputList(L, numBin-1, "Lvals.ind.mult."+Ne_str+".txt");
        }
        
        cleanUp_mt();
    }
    
    else {
        n_max = 1;
        setLengthDistribution(Tafter);
        setMaxS();
        init_mt();
        if (complete_flag){
            alpha_ind = int((max_alpha-min_alpha)/alpha_delta)+1;
            ne1_ind = int(max_ne1-min_ne1)/ne1_delta+1;
            r_ind = int((log10(max_r)-log10(min_r))/r_delta)+1;
            u_ind = int((log10(max_u)-log10(min_u))/u_delta)+1;
            //alpha_delta = (max_alpha - min_alpha)/intervals;
            //ne1_delta = (max_ne1 - min_ne1)/intervals;
            //r_delta = (log10(max_r) - log10(min_r))/intervals;
            //u_delta = (log10(max_u) - log10(min_u))/intervals;
            //tshort_delta = min(Tafter, Tbefore);
        }

        cout << "columns(index):" << endl;
        cout << "alpha,Ne1,r,u,Tbefore,Tafter,Tshort,pCorrect"<< endl;
        cout << "column intervals" << endl;
        cout << alpha_delta << "," << ne1_delta << "," << r_delta << "," << u_delta << "," << Tbefore << "," << Tafter << "," << tshort_delta << "," << endl;
        cout << "column minvals" << endl;
        cout << min_alpha << "," << min_ne1 << "," << log10(min_r) << "," << log10(min_u) << "," << min_tbefore << "," << min_tafter << "," << min_tshort << "," << endl;
        cout << "column maxvals" << endl;
        cout << max_alpha << "," << max_ne1 << "," << log10(max_r) << "," << log10(max_u) << "," << t2_max*Tbefore << "," << t1_max*Tafter << "," << max_tshort << "," << endl;
        cout << "number of col vals: " << alpha_ind << " " << ne1_ind << " " << r_ind << " " << u_ind << " " << t2_max << " " << t1_max << " " << endl;
        cout << "number of rows: " << alpha_ind*ne1_ind*r_ind*u_ind*t2_max*t1_max << endl;

        // for converting r & u with fixed precision
        stringstream ss_r, ss_u;
    
        for (int aind = 0; aind < alpha_ind; aind++){
            alpha = min_alpha + aind*alpha_delta;
            for (int neind = 0; neind < ne1_ind; neind++){
                Ne1 = min_ne1 + neind*ne1_delta;
                for (int rind = 0; rind < r_ind; rind++){
                    r = pow(10.0, (rind*r_delta)+log10(min_r));
                    for (int uind = 0; uind < u_ind; uind++){
                        u = pow(10.0, (uind*u_delta) + log10(min_u));
                        for(ct2=1; ct2 < t2_max+1; ct2++){
                            for(ct1=1; ct1 < t1_max+1; ct1++){
                                cTbefore = Tbefore*ct2;
                                cTafter_test = Tafter*ct1;
                                for (int i=0; i<2; i++){
                                    cTafter[i] = Tafter*ct1 - T_short[i]/2;
                                    cTafter_real[i] = Tafter*ct1 - T_short[i];
                                }
                                updateLengthDistribution(cTafter_test, false);
                                setGlobalVars();
                                for (int tidx=0; tidx < nthreads; tidx++){
                                    mythreads[tidx] = thread(&threadfunc_single, &thread_args[tidx]);
                                }
                                for (int tidx=0; tidx < nthreads; tidx++){
                                    mythreads[tidx].join();
                                }
                            }
                            summarize();
                            ss_r << fixed << setprecision(3) << log10(r);
                            ss_u << fixed << setprecision(3) << log10(u);
                            string indstr = to_string(aind)+","+to_string(int(Ne1))+","+ss_r.str()+","+ss_u.str()+","+
                                            to_string(int(cTbefore))+",";
                            ss_r.str(std::string());
                            ss_u.str(std::string());
                           
                            if (outpath!=""){
                                bool is_head = ( (aind==0) && (neind==0) && (rind==0) && (uind==0) && (ct2==1) );
                                outputMatrix(SumRes, t1_max, n_max, outpath, indstr, is_head);
                            }
                            else
                                printMatrix(SumRes, t1_max, n_max, false,indstr);
                            clearRes_mt();
                        }
                    }
                }
            }
        }
        cleanUp_mt();
    }
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    cout << "Runtime: " << duration.count()*1e-6 << "s" << endl;
    return EXIT_SUCCESS;
}


