#include "urqmd_facade.h"
#include <iostream>
#include <cstdlib>

//
// C++ implementation of UrQMD_main FORTRAN subroutine
// Translates the main event loop logic from FORTRAN to C++ using procedural style
//

// Define FORTRAN COMMON block name mappings
#define energies F77_NAME(energies, ENERGIES)
#define sys F77_NAME(sys, SYS)
#define rsys F77_NAME(rsys, RSYS)
#define inputs F77_NAME(inputs, INPUTS)
#define input2 F77_NAME(input2, INPUT2)
#define options F77_NAME(options, OPTIONS)
#define loptions F77_NAME(loptions, LOPTIONS)
#define stables F77_NAME(stables, STABLES)
#define comseed F77_NAME(comseed, COMSEED)
#define colltab F77_NAME(colltab, COLLTAB)
#define inewpart F77_NAME(inewpart, INEWPART)
#define coor F77_NAME(coor, COOR)
#define isys F77_NAME(isys, ISYS)
#define pots F77_NAME(pots, POTS)

// COMMON block declarations (matching FORTRAN layout)
// These allow direct access to FORTRAN COMMON blocks from C++

// /energies/ common block
extern "C" {
    struct {
        double Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau;
    } energies;
}

// /sys/ common block
extern "C" {
    struct {
        int npart, nbar, nmes, ctag, nsteps, uid_cnt;
        int ranseed, event, Ap, At, Zp, Zt, eos, dectag;
        int NHardRes, NSoftRes, NDecRes, NElColl, NBlColl;
        int success;  // logical -> int mapping
        int npartcoal, nclus;
    } sys;
}

// /rsys/ common block
extern "C" {
    struct {
        double time, acttime, bdist, bimp, bmin, ebeam, ecm;
    } rsys;
}

// /inputs/ common block
extern "C" {
    struct {
        int nevents, spityp[2], prspflg, trspflg;
        int spiso3[2], outsteps, bflag, srtflag, efuncflag, nsrt;
        int firstev, npb;
    } inputs;
}

// /input2/ common block
extern "C" {
    struct {
        double srtmin, srtmax, pbeam, betann, betatar, betapro;
        double pbmin, pbmax;
    } input2;
}

// /options/ common block
extern "C" {
    struct {
        int CTOption[400];
        double CTParam[400];
    } options;
}

// /loptions/ common block
extern "C" {
    struct {
        int fixedseed, bf13, bf14, bf15, bf16, bf19, bf20;  // logical -> int
    } loptions;
}

// /stables/ common block
extern "C" {
    struct {
        int nstable;
        int stabvec[20];
    } stables;
}

// /comseed/ common block
extern "C" {
    struct {
        int firstseed;  // logical -> int
    } comseed;
}

// /colltab/ common block
extern "C" {
    struct {
        double cttime[100001];  // 0:ncollmax
        double ctsqrts[100000];
        double ctsigtot[100000];
        double tmin;
        int cti1[100000];
        int cti2[100000];
        int nct;
        int actcol;
        int ctvalid[100000];  // logical -> int
        int ctsav[100000];
        int nsav;
        int apt;
        double ctcolfluc[100000];
    } colltab;
}

// /inewpart/ common block
extern "C" {
    struct {
        int itypnew[1000];
        int i3new[1000];
        int itot[1000];
        int inew[1000];
        int nexit;
        int iline;
        int pslot[2];
        int nstring1, nstring2;
        int itypold[2];
        int iso3old[2];
    } inewpart;
}

// /coor/ common block (particle coordinates)
extern "C" {
    struct {
        double r0[100000];
        double rx[100000];
        double ry[100000];
        double rz[100000];
        double p0[100000];
        double px[100000];
        double py[100000];
        double pz[100000];
        double fmass[100000];
        double rww[100000];
        double dectime[100000];
    } coor;
}

// /isys/ common block (particle properties)
extern "C" {
    struct {
        int spin[100000];
        int ncoll[100000];
        int charge[100000];
        int ityp[100000];
        int lstcoll[100000];
        int iso3[100000];
        int origin[100000];
        int uid[100000];
    } isys;
}

// /pots/ common block
extern "C" {
    struct {
        double Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky;
        double gamYuk, drPau, dpPau, gw, sgw, delr, fdel;
        double dt, da, db, dtimestep;
    } pots;
}

// Parameters
const double emnuc = 0.938;  // nucleon mass

// Macro to access COMMON blocks (using the defined names)
#define SYS sys
#define RSYS rsys
#define INPUTS inputs
#define OPTIONS options
#define LOPTIONS loptions
#define COMSEED comseed
#define COLLTAB colltab
#define INEWPART inewpart
#define COOR coor
#define ISYS isys
#define POTS pots
#define STABLES stables

//
// Main UrQMD function - C++ translation of UrQMD_main subroutine
//
void UrQMD_main()
{
    // Local variables matching FORTRAN declarations
    int i, j, k, steps, ii, ocharge, ncharge, mc, mp, noc, it1, it2;
    double sqrts_val, otime, xdummy, st;
    bool isstable;
    int stidx, CTOsave;
    int cti1sav, cti2sav;
    
    // Hydro variables
    double thydro_start, thydro, nucrad_val;
    bool lhydro;
    
    //
    // Numerical/technical initialization
    //
    int init_flag = 0;
    uinit(&init_flag);
    
    //
    // Main program
    //
    mc = 0;
    mp = 0;
    noc = 0;
    
    //
    // Loop over all events
    //
    for (SYS.event = 1; SYS.event <= INPUTS.nevents; SYS.event++) {
        
        // Start event here
        
        // Declare variables that might be used with goto statements
        int one_val = 1;
        int zero_val = 0;
        int neg_one_val = -1;
        
        // Time is the system time at the BEGINNING of every timestep
        RSYS.time = 0.0;
        
        // Hydro flag, hydro should be called only once
        lhydro = true;
        
        // Initialize random number generator
        // Call auto-seed generator only for first event and if no seed was fixed
        if (!COMSEED.firstseed && !LOPTIONS.fixedseed) {
            SYS.ranseed = -(1 * abs(SYS.ranseed));
            sseed(&SYS.ranseed);
        } else {
            COMSEED.firstseed = 0;  // false
        }
        
        printf("event# %d %d\n", SYS.event, SYS.ranseed);
        
        // Init resonance reconstruction (only f13)
        if (!LOPTIONS.bf13 && OPTIONS.CTOption[68-1] == 1) {
            init_resrec();
        }
        
        //
        // Initialization of physics quantities
        //
        init();
        init_eccentricity();
        
        // If we are reading old events, check the success of the read-in
        if (OPTIONS.CTOption[40-1] != 0 && !SYS.success) {
            break;  // exit event loop
        }
        
        // Hydro switch
        if (OPTIONS.CTOption[45-1] == 1) {
            // Hydro start time (nuclei have passed each other)
            // ebeam is only the kinetic energy
            // CTParam(65) is useful for the variation of the start time
            // default value is one
            nucrad_val = nucrad(&SYS.Ap);
            thydro_start = OPTIONS.CTParam[65-1] * 2.0 * nucrad_val * 
                          sqrt(2.0 * emnuc / RSYS.ebeam);
            printf("hydro starts after %.2f fm/c\n", thydro_start);
            
            // Lower limit for hydro start time
            if (thydro_start < OPTIONS.CTParam[63-1]) {
                thydro_start = OPTIONS.CTParam[63-1];
                printf("... extended to %.2f fm/c\n", OPTIONS.CTParam[63-1]);
            }
        }
        
        // Old time if an old fort.14 is used
        if (OPTIONS.CTOption[40-1] != 0) {
            RSYS.time = RSYS.acttime;
        }
        
        // Output preparation
        
        // Write headers to file
        int unit;
        unit = 13; output(&unit);
        unit = 14; output(&unit);
        unit = 15; output(&unit);
        unit = 16; output(&unit);
        
        if (SYS.event == 1) {
            unit = 17; output(&unit);
            osc_header();
            osc99_header();
        }
        
        osc99_event(&neg_one_val);
        
        // For CTOption(4)=1 : output of initialization configuration
        if (OPTIONS.CTOption[4-1] == 1) {
            file14out(&zero_val);
        }
        
        // Participant/spectator model
        if (OPTIONS.CTOption[28-1] != 0) {
            double x1 = 0.5 * RSYS.bimp;
            double x2 = -(0.5 * RSYS.bimp);
            rmspec(&x1, &x2);
        }
        
        // Compute time of output
        otime = INPUTS.outsteps * POTS.dtimestep;
        
        // Reset time step counter
        steps = 0;
        
        //
        // Loop over all timesteps
        //
        for (steps = 1; steps <= SYS.nsteps; steps++) {
            
            // Declare all variables that might be used with goto statements
            int idx1, idx2;
            double sqrts_val_local;
            
            // Store coordinates in arrays with *_t
            // This is needed for MD type propagation
            // NOTE: MD propagation requires r0_t, rx_t, ry_t, rz_t temporary arrays
            // which are in the /mdprop/ COMMON block (see coms.inc)
            // This functionality is preserved but not implemented in this minimal facade
            if (SYS.eos != 0) {
                for (j = 0; j < SYS.npart; j++) {
                    // TODO: Implement MD coordinate storage when needed
                    // r0_t(j) = r0(j); rx_t(j) = rx(j); etc.
                }
            }
            
            // We are at the beginning of the timestep, set current time (acttime)
            RSYS.acttime = RSYS.time;
            
            // Option for MD without collision term
            if (OPTIONS.CTOption[16-1] != 0) {
                goto label_103;
            }
            
            // Load collision table with next collisions in current timestep
            colload();
            
            // Check for collisions in time-step, nct = # of collisions in table
            if (COLLTAB.nct > 0) {
                
                // Entry-point for collision loop in case of full colload after every coll.
            label_101:
                k = 0;
                
                // Normal entry-point for collision loop
            label_100:
                
                // Get next collision
                getnext(&k);
                
                // Exit collision loop if no collisions are left
                if (k == 0) goto label_102;
                
                // hp call hydro if start time is reached
                if (OPTIONS.CTOption[45-1] == 1) {
                    if (COLLTAB.cttime[k] > thydro_start && lhydro) {
                        
                        if (OPTIONS.CTOption[62-1] == 1) {
                            prepout();
                            int zero = 0;
                            file14out(&zero);
                            restore();
                        }
                        
                        st = thydro_start - RSYS.acttime;
                        cascstep(&RSYS.acttime, &st);
                        
                        // hp all particle arrays will be modified by hydro
                        printf("starting hydro\n");
                        hydro(&thydro_start, &thydro);
                        RSYS.acttime = thydro_start;
                        lhydro = false;
                        
                        if (OPTIONS.CTOption[50-1] == 1) goto label_10_continue;
                        
                        if (thydro > 1.0e-8 || OPTIONS.CTOption[48-1] == 1) {
                            // hp full update of collision table
                            colload();
                            goto label_101;
                        }
                    }
                }
                
                // Propagate all particles to next collision time
                // Store actual time in acttime, propagation time st=cttime(k)-acttime
                st = COLLTAB.cttime[k] - RSYS.acttime;
                cascstep(&RSYS.acttime, &st);
                
                // New actual time (for upcoming collision)
                RSYS.acttime = COLLTAB.cttime[k];
                
                // Perform collision
                idx1 = COLLTAB.cti1[k-1];
                idx2 = COLLTAB.cti2[k-1];
                
                if (idx2 > 0) {
                    sqrts_val_local = sqrts(&idx1, &idx2);
                    if (fabs(sqrts_val_local - COLLTAB.ctsqrts[k-1]) > 1.0e-3) {
                        fprintf(stderr, " ***(E) wrong collision update (col) ***\n");
                        fprintf(stderr, "%d %d %d %.6f %.6f\n", k, idx1, idx2,
                               COLLTAB.ctsqrts[k-1], sqrts_val_local);
                    }
                } else if (idx2 == 0) {
                    if (fabs(COOR.fmass[idx1-1] - COLLTAB.ctsqrts[k-1]) > 1.0e-3) {
                        fprintf(stderr, " *** main(W) wrong collision update (decay)\n");
                        fprintf(stderr, "%d %d %d %.6e %.6f %.6f\n", 
                               SYS.ctag, idx1, ISYS.ityp[idx1-1],
                               COOR.dectime[idx1-1], COOR.fmass[idx1-1],
                               COLLTAB.ctsqrts[k-1]);
                    }
                }
                
                ocharge = ISYS.charge[idx1-1];
                if (idx2 > 0) ocharge += ISYS.charge[idx2-1];
                
                // Store quantities in local variables for charge conservation check
                it1 = ISYS.ityp[idx1-1];
                if (idx2 > 0) it2 = ISYS.ityp[idx2-1];
                
                // Increment "dirty" collision counter
                if (idx2 > 0) {  // scatter
                    mc++;
                }
                
                // Perform scattering/decay
                cti1sav = COLLTAB.cti1[k-1];
                cti2sav = COLLTAB.cti2[k-1];
                
                scatter(&COLLTAB.cti1[k-1], &COLLTAB.cti2[k-1],
                                          &COLLTAB.ctsigtot[k-1], &COLLTAB.ctsqrts[k-1],
                                          &COLLTAB.ctcolfluc[k-1]);
                
                //
                // Update collision table
                //
                // Normal update mode
                if (OPTIONS.CTOption[17-1] == 0) {
                    if (INEWPART.nexit == 0) {
                        // New collision partners for pauli-blocked states (nexit=0)
                        if (COLLTAB.cti1[k-1] != cti1sav || COLLTAB.cti2[k-1] != cti2sav) {
                            COLLTAB.cti1[k-1] = cti1sav;
                            COLLTAB.cti2[k-1] = cti2sav;
                        }
                        int one = 1;
                        collupd(&COLLTAB.cti1[k-1], &one);
                        if (COLLTAB.cti2[k-1] > 0) {
                            collupd(&COLLTAB.cti2[k-1], &one);
                        }
                    } else {
                        ncharge = 0;
                        // New collision partners for scattered/produced particles (nexit><0)
                        for (i = 0; i < INEWPART.nexit; i++) {
                            // ncharge is used for charge conservation check
                            ncharge += ISYS.charge[INEWPART.inew[i]-1];
                            int one = 1;
                            collupd(&INEWPART.inew[i], &one);
                        }
                        
                        // Charge conservation check
                        if (ocharge != ncharge) {
                            fprintf(stderr, "ch-conservation error coll/dec %d\n", SYS.ctag);
                            fprintf(stderr, "   it1: %d   it2: %d\n", it1, it2);
                            fprintf(stderr, "   ch: %d %d\n", ocharge, ncharge);
                            fprintf(stderr, "cti1(k),cti2(k),ctsigtot(k),ctsqrts(k)\n");
                            fprintf(stderr, "%d %d %.6f %.6f\n", 
                                   COLLTAB.cti1[k-1], COLLTAB.cti2[k-1],
                                   COLLTAB.ctsigtot[k-1], COLLTAB.ctsqrts[k-1]);
                        }
                    }
                    
                    // Update collisions for partners of annihilated particles
                    for (ii = 0; ii < COLLTAB.nsav; ii++) {
                        int one = 1;
                        collupd(&COLLTAB.ctsav[ii], &one);
                    }
                    COLLTAB.nsav = 0;
                } else {  // (CTOption(17).ne.0)
                    // Full collision load
                    colload();
                }
                
                if (OPTIONS.CTOption[17-1] == 0) goto label_100;
                goto label_101;
                
                // This is the point to jump to after all collisions in the timestep
                // have been taken care of
            label_102:
                ;  // continue
                
            }  // (nct.gt.0)
            
            // After all collisions in the timestep are done, propagate to end of
            // the timestep.
            
            // Point to jump to in case of MD without collision term
        label_103:
            
            // Increment timestep
            RSYS.time = RSYS.time + POTS.dtimestep;
            
            // After all collisions in the timestep are done, propagate to end of
            // the timestep.
            st = RSYS.time - RSYS.acttime;
            cascstep(&RSYS.acttime, &st);
            
            // In case of potential interaction, do MD propagation step
            if (SYS.eos != 0) {
                // Set initial conditions for MD propagation-step
                // NOTE: This requires the /mdprop/ COMMON block with temporary arrays
                // r0_t, rx_t, ry_t, rz_t (see coms.inc line 96-98)
                for (j = 0; j < SYS.npart; j++) {
                    // TODO: Implement MD coordinate restoration when needed
                    // r0(j) = r0_t(j); rx(j) = rx_t(j); etc.
                }
                
                // Now molecular dynamics trajectories
                proprk(&RSYS.time, &POTS.dtimestep);
            }
            
            // Perform output if desired
            if ((steps % INPUTS.outsteps) == 0 && steps < SYS.nsteps) {
                if (OPTIONS.CTOption[28-1] == 2) {
                    spectrans(&otime);
                }
                if (OPTIONS.CTOption[62-1] == 1) {
                    prepout();
                }
                file14out(&steps);
                if (OPTIONS.CTOption[64-1] == 1) {
                    file13out(&steps);
                }
                if (OPTIONS.CTOption[62-1] == 1) {
                    restore();
                    colload();
                }
                if (OPTIONS.CTOption[55-1] == 1) {
                    osc_vis(&steps);
                }
            }  // output handling
            
        }  // time step loop
        
        // e
        RSYS.acttime = RSYS.time;
        
        // Optional decay of all unstable particles before final output
        // DANGER: pauli-blocked decays are not performed !!!
        if (OPTIONS.CTOption[18-1] == 0) {
            // No do-loop is used because npart changes in loop-structure
            i = 0;
            COLLTAB.nct = 0;
            COLLTAB.actcol = 0;
            
            // Disable Pauli-Blocker for final decays
            CTOsave = OPTIONS.CTOption[10-1];
            OPTIONS.CTOption[10-1] = 1;
            
            // Decay loop structure starts here
        label_40:
            i++;
            
            // Is particle unstable
            if (COOR.dectime[i-1] < 1.0e30) {
            label_41:
                isstable = false;
                for (stidx = 0; stidx < STABLES.nstable; stidx++) {
                    if (ISYS.ityp[i-1] == STABLES.stabvec[stidx]) {
                        isstable = true;
                    }
                }
                if (!isstable) {
                    // Perform decay
                    int zero = 0;
                    double zero_d = 0.0;
                    scatter(&i, &zero, &zero_d, 
                                              &COOR.fmass[i-1], &xdummy);
                    
                    // Backtracing if decay-product is unstable itself
                    if (COOR.dectime[i-1] < 1.0e30) goto label_41;
                }
            }
            
            // Check next particle
            if (i < SYS.npart) goto label_40;
            
        }  // final decay
        
        OPTIONS.CTOption[10-1] = CTOsave;
        
        // Final output
        if (OPTIONS.CTOption[64-1] == 1) {
            coalescence();
        }
        
        file13out(&SYS.nsteps);
        if (OPTIONS.CTOption[50-1] == 0) {
            file14out(&SYS.nsteps);
        }
        file16out();
        if (OPTIONS.CTOption[50-1] == 0 && OPTIONS.CTOption[55-1] == 0) {
            osc_event();
        }
        if (OPTIONS.CTOption[50-1] == 0 && OPTIONS.CTOption[55-1] == 1) {
            osc_vis(&SYS.nsteps);
        }
        osc99_event(&one_val);
        osc99_eoe();
        
        mp += SYS.npart;
        if (SYS.ctag == 0) {
            printf("(W) No collision in event %d\n", SYS.event);
            noc++;
        }
        
        // End of event loop
    label_10_continue:
        continue;
    }  // event loop
    
    printf("no. of collisions = %.1f (per event)\n", 
           mc / static_cast<double>(INPUTS.nevents));
    printf("final particles   = %.1f (per event)\n",
           mp / static_cast<double>(INPUTS.nevents));
    printf("empty events      : %d = %.1f%%\n",
           noc, noc * 100.0 / static_cast<double>(INPUTS.nevents));
}
