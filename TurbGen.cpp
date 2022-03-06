/***********************************************
 *** Generate turbulent field with TurbGen.h ***
 *** Written by Christoph Federrath, 2022.   ***
************************************************/

#include <string>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include "TurbGen.h"

using namespace std;

// global variables
enum {X, Y, Z};
static const string ProgSign = "TurbGen: ";
int verbose = 1; // 0: all standard output off (quiet mode)
int ndim = 3; // dimensionality (must be 1 or 2 or 3)
int N[3] = {64, 64, 64}; // number of cells in turbulent output field in x, y, z
double L[3] = {1.0, 1.0, 1.0}; // size of box in x, y, z
double k_min = 2.0; // minimum wavenumber for turbulent field (in units of 2pi / L[X])
double k_max = 20.0; // minimum wavenumber for turbulent field (in units of 2pi / L[X])
int spect_form = 2; // 0: band/rectangle/constant, 1: paraboloid, 2: power law
double power_law_exp = -2.0; // if spect_form == 2: power-law exponent (e.g., -2 would be Burgers,
                             // or -5/3 would be Kolmogorov, or 1.5 would Kazantsev)
double angles_exp = 1.0; // if spect_form == 2: spectral sampling of angles;
                         // number of modes (angles) in k-shell surface increases as k^angles_exp.
                         // For full sampling, angles_exp = 2.0; for healpix-type sampling, angles_exp = 0.0.
double sol_weight = 0.5; // solenoidal weight: 1.0: solenoidal driving, 0.0: compressive driving, 0.5: natural mixture
int random_seed = 140281; // random seed for this turbulent realisation


// forward functions
int ParseInputs(const vector<string> Argument);
void HelpMe(void);


int main(int argc, char * argv[])
{
    /// Parse inputs
    vector<string> Arguments(argc);
    for (int i = 0; i < argc; i++) Arguments[i] = static_cast<string>(argv[i]);
    if (ParseInputs(Arguments) == -1)
    {
        if (verbose > 0) cout<<endl<<"Error in ParseInputs(). Exiting."<<endl;
        HelpMe();
        return -1;
    }

    if (verbose > 0) cout<<ProgSign+"started..."<<endl;

    long starttime = time(NULL);

    // create TurbGen class object
    TurbGen tg = TurbGen();
    tg.set_verbose(verbose);

    // set dimensionality dependencies
    if (ndim < 3) { N[Z] = 1; L[Z] = 1.0; } // 2D
    if (ndim < 2) { N[Y] = 1; L[Y] = 1.0; } // 1D

    // define cell size of uniform grid
    double del[3];
    for (int d = 0; d < 3; d++) del[d] = L[d] / N[d];
    // always generate within [0,L], cell-centered; user can shift output to target physical location, if needed
    double pos_beg[3], pos_end[3]; // start and end coordinates of output grid
    for (int d = 0; d < 3; d++) {
        pos_beg[d] = del[d]/2.0; // first cell coordinate (cell center)
        pos_end[d] = L[d] - del[d]/2.0; // last cell coordinate (cell center)
    }
    // total number of output grid cells
    long ntot = 1; for (int d = 0; d < 3; d++) ntot *= N[d];
    // output grid, which receives the turbulent field (up to ndim = 3)
    float * grid_out[3];
    // allocate
    for (int d = 0; d < ndim; d++) grid_out[d] = new float[ntot];

    // initialise generator to return a single turbulent realisation based on input parameters
    tg.init_single_realisation(ndim, L, k_min, k_max, spect_form, power_law_exp, angles_exp, sol_weight, random_seed);

    // call to return uniform grid(s) with ndim components of the turbulent field at requested positions
    tg.get_turb_vector_unigrid(pos_beg, pos_end, N, grid_out);

    // indices in x,y,z, i.e., i,j,k, and 1D index
    int i, j, k; long index;
    i = 3; j = 0; k = 10;
    index = k*N[Y]*N[X] + j*N[X] + i;
    printf("e.g., x-component of turbulent vector field at index i,j,k = (%i %i %i) is %f\n", i,j,k, grid_out[0][index]);
    i = 9; j = 1; k = 4;
    index = k*N[Y]*N[X] + j*N[X] + i;
    printf("e.g., z-component of turbulent vector field at index i,j,k = (%i %i %i) is %f\n", i,j,k, grid_out[2][index]);

    // compute mean and RMS of generated turbulent field and then re-normalise to mean=0 and std=1
    double mean[3] = {0.0, 0.0, 0.0};
    double rms[3] = {0.0, 0.0, 0.0};
    double std[3] = {0.0, 0.0, 0.0};
    for (int d = 0; d < ndim; d++) {
        for (long ni = 0; ni < ntot; ni++) {
            mean[d] += grid_out[d][ni];
            rms [d] += pow(grid_out[d][ni],2.0);
        }
        mean[d] /= ntot; // mean
        rms[d] = sqrt(rms[d] / ntot); // root mean squared
        std[d] = sqrt(rms[d]*rms[d] - mean[d]*mean[d]); // standard deviation
    }
    // re-normalise to mean=0 and std=1
    for (int d = 0; d < ndim; d++) {
        for (long ni = 0; ni < ntot; ni++) {
            grid_out[d][ni] -= mean[d];
            grid_out[d][ni] /= std[d];
        }
    }
    // re-compute mean and std
    for (int d = 0; d < ndim; d++) {
        mean[d] = 0.0; rms[d] = 0.0;
        for (long ni = 0; ni < ntot; ni++) {
            mean[d] += grid_out[d][ni];
            rms [d] += pow(grid_out[d][ni],2.0);
        }
        mean[d] /= ntot; // mean
        rms[d] = sqrt(rms[d] / ntot); // root mean squared
        std[d] = sqrt(rms[d]*rms[d] - mean[d]*mean[d]); // standard deviation
    }
    printf("Turbulent vector field mean (x,y,z) = (%e %e %e)\n", mean[0], mean[1], mean[2]);
    printf("Turbulent vector field  rms (x,y,z) = (%e %e %e)\n",  rms[0],  rms[1],  rms[2]);
    printf("Turbulent vector field  std (x,y,z) = (%e %e %e)\n",  std[0],  std[1],  std[2]);
    printf("Turbulent vector field total 3D std = %e\n", sqrt(std[0]*std[0]+std[1]*std[1]+std[2]*std[2]));

    // clean up
    for (int d = 0; d < ndim; d++) {
      if (grid_out[d]) delete [] grid_out[d]; grid_out[d] = NULL;
    }

    /// print out wallclock time used
    long endtime = time(NULL);
    long duration = endtime-starttime;
    if (verbose > 0) cout<<ProgSign+"Total runtime: "<<duration<<"s"<<endl;

  return 0;
}



/** ------------------------- ParseInputs ----------------------------
 **  Parses the command line Arguments
 ** ------------------------------------------------------------------ */
int ParseInputs(const vector<string> Argument)
{
    static const string FuncSign = ProgSign+"ParseInputs: ";
    stringstream dummystream;

    /// read tool specific options
    if (Argument.size() < 2)
    {
        cout << endl << FuncSign+"Specify at least 1 argument." << endl;
        return -1;
    }

    for (unsigned int i = 1; i < Argument.size(); i++)
    {
        if (Argument[i] != "" && Argument[i] == "-h")
        {
            HelpMe(); exit(0);
        }
        if (Argument[i] != "" && Argument[i] == "-verbose")
        {
            if (Argument.size()>i+1) {
                dummystream << Argument[i+1]; dummystream >> verbose; dummystream.clear();
            } else return -1;
        }
        if (Argument[i] != "" && Argument[i] == "-ndim")
        {
            if (Argument.size()>i+1) {
                dummystream << Argument[i+1]; dummystream >> ndim; dummystream.clear();
            } else return -1;
        }
        if (Argument[i] != "" && Argument[i] == "-N")
        {
            if (Argument.size()>i+ndim) {
                dummystream << Argument[i+1]; dummystream >> N[X]; dummystream.clear();
                if (ndim > 1) { dummystream << Argument[i+2]; dummystream >> N[Y]; dummystream.clear(); }
                if (ndim > 2) { dummystream << Argument[i+3]; dummystream >> N[Z]; dummystream.clear(); }
            } else return -1;
        }
        if (Argument[i] != "" && Argument[i] == "-L")
        {
            if (Argument.size()>i+ndim) {
                dummystream << Argument[i+1]; dummystream >> L[X]; dummystream.clear();
                if (ndim > 1) { dummystream << Argument[i+2]; dummystream >> L[Y]; dummystream.clear(); }
                if (ndim > 2) { dummystream << Argument[i+3]; dummystream >> L[Z]; dummystream.clear(); }
            } else return -1;
        }
        if (Argument[i] != "" && Argument[i] == "-kmin")
        {
            if (Argument.size()>i+1) {
                dummystream << Argument[i+1]; dummystream >> k_min; dummystream.clear();
            } else return -1;
        }
        if (Argument[i] != "" && Argument[i] == "-kmax")
        {
            if (Argument.size()>i+1) {
                dummystream << Argument[i+1]; dummystream >> k_max; dummystream.clear();
            } else return -1;
        }
        if (Argument[i] != "" && Argument[i] == "-spect_form")
        {
            if (Argument.size()>i+1) {
                dummystream << Argument[i+1]; dummystream >> spect_form; dummystream.clear();
            } else return -1;
        }
        if (Argument[i] != "" && Argument[i] == "-power_law_exp")
        {
            if (Argument.size()>i+1) {
                dummystream << Argument[i+1]; dummystream >> power_law_exp; dummystream.clear();
            } else return -1;
        }
        if (Argument[i] != "" && Argument[i] == "-angles_exp")
        {
            if (Argument.size()>i+1) {
                dummystream << Argument[i+1]; dummystream >> angles_exp; dummystream.clear();
            } else return -1;
        }
        if (Argument[i] != "" && Argument[i] == "-sol_weight")
        {
            if (Argument.size()>i+1) {
                dummystream << Argument[i+1]; dummystream >> sol_weight; dummystream.clear();
            } else return -1;
        }
        if (Argument[i] != "" && Argument[i] == "-random_seed")
        {
            if (Argument.size()>i+1) {
                dummystream << Argument[i+1]; dummystream >> random_seed; dummystream.clear();
            } else return -1;
        }

    } // loop over all args

    /// print out parsed values
    if (verbose > 1) {
        cout << FuncSign+"Command line arguments: ";
        for (unsigned int i = 0; i < Argument.size(); i++) cout << Argument[i] << " ";
        cout << endl;
    }
    return 0;

} // end: ParseInputs()


/** --------------------------- HelpMe -------------------------------
 **  Prints out usage information
 ** ------------------------------------------------------------------ */
void HelpMe(void)
{
    if (verbose > 0) {
        cout << endl
        << ProgSign+"Generates a turbulent field of mean=0 and std=1, with specified parameters (see OPTIONS)." << endl << endl
        << "Syntax:" << endl
        << " TurbGen [<OPTIONS>]" << endl << endl
        << "   <OPTIONS>:" << endl
        << "     -ndim <1, 2, 3>        : number of spatial dimensions; (default: 3)" << endl
        << "     -N <Nx [Ny [Nz]]>      : number of grid cells for generated turbulent field; (default: 64 64 64)" << endl
        << "     -L <Lx [Ly [Lz]]>      : physical size of box in which to generate turbulent field; (default: 1.0 1.0 1.0)" << endl
        << "     -kmin <val>            : minimum wavenumber of generated field in units of 2pi/L[X]; (default: 2.0)" << endl
        << "     -kmax <val>            : maximum wavenumber of generated field in units of 2pi/L[X]; (default: 20.0)" << endl
        << "     -spect_form <0, 1, 2>  : spectral form: 0 (band/rectangle/constant), 1 (paraboloid), 2 (power law); (default: 2)" << endl
        << "     -power_law_exp <val>   : if spect_form 2: power-law exponent of power-law spectrum (e.g. Kolmogorov: -5/3, Burgers: -2, Kazantsev: 1.5); (default -2.0)" << endl
        << "     -angles_exp <val>      : if spect_form 2: angles exponent for sparse sampling (e.g., 2.0: full sampling, 0.0: healpix-like sampling); (default 1.0)" << endl
        << "     -sol_weight <val>      : solenoidal weight: 1.0 (divergence-free field), 0.5 (natural mix), 0.0 (curl-free field); (default 0.5)" << endl
        << "     -random_seed <val>     : random seed for turbulent field; (default 140281)" << endl
        << "     -verbose <0, 1, 2>     : 0 (no shell output), 1 (standard shell output), 2 (more shell output); (default 1)" << endl
        << "     -h                     : print this help message" << endl
        << endl
        << "Example: TurbGen -ndim 2 -L 1.0 1.0"
        << endl << endl;
    }
}
