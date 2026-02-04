#ifndef __urqmd_facade_h__
#define __urqmd_facade_h__

#include "StarCallf77.h"
#include <cstdio>
#include <cmath>

//
// C++ facade for UrQMD FORTRAN subroutine UrQMD_main
// This provides a procedural C++ interface to the FORTRAN UrQMD event generator
//

// Define FORTRAN name mappings using F77_NAME macro
#define uinit F77_NAME(uinit, UINIT)
#define sseed F77_NAME(sseed, SSEED)
#define init F77_NAME(init, INIT)
#define init_resrec F77_NAME(init_resrec, INIT_RESREC)
#define init_eccentricity F77_NAME(init_eccentricity, INIT_ECCENTRICITY)

#define output F77_NAME(output, OUTPUT)
#define osc_header F77_NAME(osc_header, OSC_HEADER)
#define osc99_header F77_NAME(osc99_header, OSC99_HEADER)
#define osc99_event F77_NAME(osc99_event, OSC99_EVENT)
#define osc99_eoe F77_NAME(osc99_eoe, OSC99_EOE)
#define osc_event F77_NAME(osc_event, OSC_EVENT)
#define osc_vis F77_NAME(osc_vis, OSC_VIS)
#define file13out F77_NAME(file13out, FILE13OUT)
#define file14out F77_NAME(file14out, FILE14OUT)
#define file16out F77_NAME(file16out, FILE16OUT)

#define colload F77_NAME(colload, COLLOAD)
#define getnext F77_NAME(getnext, GETNEXT)
#define cascstep F77_NAME(cascstep, CASCSTEP)
#define scatter F77_NAME(scatter, SCATTER)
#define collupd F77_NAME(collupd, COLLUPD)
#define proprk F77_NAME(proprk, PROPRK)

#define rmspec F77_NAME(rmspec, RMSPEC)
#define prepout F77_NAME(prepout, PREPOUT)
#define restore F77_NAME(restore, RESTORE)
#define hydro F77_NAME(hydro, HYDRO)
#define spectrans F77_NAME(spectrans, SPECTRANS)
#define coalescence F77_NAME(coalescence, COALESCENCE)

#define nucrad F77_NAME(nucrad, NUCRAD)
#define sqrts F77_NAME(sqrts, SQRTS)

// Forward declarations for FORTRAN subroutines
extern "C" {
    // Initialization and utility functions
    void type_of_call uinit(int *flag);
    void type_of_call sseed(int *seed);
    void type_of_call init();
    void type_of_call init_resrec();
    void type_of_call init_eccentricity();
    
    // Output functions
    void type_of_call output(int *unit);
    void type_of_call osc_header();
    void type_of_call osc99_header();
    void type_of_call osc99_event(int *flag);
    void type_of_call osc99_eoe();
    void type_of_call osc_event();
    void type_of_call osc_vis(int *steps);
    void type_of_call file13out(int *steps);
    void type_of_call file14out(int *steps);
    void type_of_call file16out();
    
    // Collision and propagation functions
    void type_of_call colload();
    void type_of_call getnext(int *k);
    void type_of_call cascstep(double *acttime, double *st);
    void type_of_call scatter(int *i1, int *i2, double *sigtot, double *sqrts, double *colfluc);
    void type_of_call collupd(int *index, int *flag);
    void type_of_call proprk(double *time, double *dtimestep);
    
    // Special functions
    void type_of_call rmspec(double *x1, double *x2);
    void type_of_call prepout();
    void type_of_call restore();
    void type_of_call hydro(double *tstart, double *tend);
    void type_of_call spectrans(double *otime);
    void type_of_call coalescence();
    
    // Utility functions
    double type_of_call nucrad(int *A);
    double type_of_call sqrts(int *i1, int *i2);
}

// Main UrQMD entry point - C++ implementation of UrQMD_main subroutine
void UrQMD_main();

#endif
