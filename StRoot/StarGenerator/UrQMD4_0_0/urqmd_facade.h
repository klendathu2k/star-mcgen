#ifndef __urqmd_facade_h__
#define __urqmd_facade_h__

#include "StarCallf77.h"
#include <cstdio>
#include <cmath>

//
// C++ facade for UrQMD FORTRAN subroutine UrQMD_main
// This provides a procedural C++ interface to the FORTRAN UrQMD event generator
//

// Forward declarations for FORTRAN subroutines that will be called
extern "C" {
    // Initialization and utility functions
    void F77_NAME(uinit, UINIT)(int *flag);
    void F77_NAME(sseed, SSEED)(int *seed);
    void F77_NAME(init, INIT)();
    void F77_NAME(init_resrec, INIT_RESREC)();
    void F77_NAME(init_eccentricity, INIT_ECCENTRICITY)();
    
    // Output functions
    void F77_NAME(output, OUTPUT)(int *unit);
    void F77_NAME(osc_header, OSC_HEADER)();
    void F77_NAME(osc99_header, OSC99_HEADER)();
    void F77_NAME(osc99_event, OSC99_EVENT)(int *flag);
    void F77_NAME(osc99_eoe, OSC99_EOE)();
    void F77_NAME(osc_event, OSC_EVENT)();
    void F77_NAME(osc_vis, OSC_VIS)(int *steps);
    void F77_NAME(file13out, FILE13OUT)(int *steps);
    void F77_NAME(file14out, FILE14OUT)(int *steps);
    void F77_NAME(file16out, FILE16OUT)();
    
    // Collision and propagation functions
    void F77_NAME(colload, COLLOAD)();
    void F77_NAME(getnext, GETNEXT)(int *k);
    void F77_NAME(cascstep, CASCSTEP)(double *acttime, double *st);
    void F77_NAME(scatter, SCATTER)(int *i1, int *i2, double *sigtot, double *sqrts, double *colfluc);
    void F77_NAME(collupd, COLLUPD)(int *index, int *flag);
    void F77_NAME(proprk, PROPRK)(double *time, double *dtimestep);
    
    // Special functions
    void F77_NAME(rmspec, RMSPEC)(double *x1, double *x2);
    void F77_NAME(prepout, PREPOUT)();
    void F77_NAME(restore, RESTORE)();
    void F77_NAME(hydro, HYDRO)(double *tstart, double *tend);
    void F77_NAME(spectrans, SPECTRANS)(double *otime);
    void F77_NAME(coalescence, COALESCENCE)();
    
    // Utility functions
    double F77_NAME(nucrad, NUCRAD)(int *A);
    double F77_NAME(sqrts, SQRTS)(int *i1, int *i2);
}

// Main UrQMD entry point - C++ implementation of UrQMD_main subroutine
void UrQMD_main();

#endif
