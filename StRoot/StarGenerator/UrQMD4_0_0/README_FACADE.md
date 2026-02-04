# UrQMD C++ Facade

This directory contains a C++ procedural facade for the FORTRAN UrQMD event generator.

## Files Created

### urqmd_facade.h
Header file containing:
- Forward declarations for all FORTRAN subroutines called by UrQMD_main
- Function prototype for the main `UrQMD_main()` entry point
- Include of StarCallf77.h for FORTRAN name mangling macros

### urqmd_facade.cxx
Implementation file containing:
- C++ translation of the FORTRAN `UrQMD_main` subroutine from urqmd.F
- COMMON block declarations as extern C structs matching FORTRAN memory layout
- Procedural implementation preserving the original FORTRAN structure and logic
- Main event loop with collision handling, propagation, and output

### address.F
FORTRAN file providing C-compatible access to COMMON blocks:
- Functions returning pointers to FORTRAN COMMON block addresses
- Updated from UrQMD3_3_1 version with changes for UrQMD 4.0.0:
  - nmax increased from 40000 to 100000 particles
  - Added success, npartcoal, nclus to /sys/ common block
  - Removed strid from /isys/ common block
- Uses iso_c_binding for C interoperability

## Structure and Design

The C++ facade follows a **procedural programming style** without heavy object-oriented encapsulation:

1. **Direct COMMON block access**: Uses extern C struct declarations to directly access FORTRAN COMMON blocks
2. **Procedural function calls**: Calls FORTRAN subroutines directly through C linkage
3. **Structure preservation**: Maintains the same control flow and logic as the original FORTRAN code
4. **Goto translation**: Converts FORTRAN goto statements to C++ labels while respecting C++ scoping rules

## Key Differences from FORTRAN

### Variable Declarations
- Variables used across goto statements are declared at the top of their scope
- Local variables (e.g., `one_val`, `zero_val`) replace inline literal declarations

### Array Indexing
- FORTRAN 1-based indexing converted to C++ 0-based indexing where accessing arrays
- COMMON block array accesses use appropriate index adjustments (e.g., `[i-1]` for element `i`)

### Logical Types
- FORTRAN LOGICAL types mapped to C++ int (0=false, non-zero=true)

### COMMON Blocks
The facade accesses FORTRAN COMMON blocks through extern C declarations:
- /sys/ - System parameters (npart, event, etc.)
- /rsys/ - Real system parameters (time, energy, etc.)
- /options/ - Control options and parameters
- /inputs/ - Input parameters
- /colltab/ - Collision table
- /coor/ - Particle coordinates
- /isys/ - Particle properties
- /pots/ - Potential parameters
- And others as needed

## Usage

To use the C++ facade:

```cpp
#include "urqmd_facade.h"

// Initialize UrQMD (if needed, via other initialization functions)
// ...

// Run the main event loop
UrQMD_main();
```

The facade integrates with the existing UrQMD FORTRAN library and requires:
- The UrQMD FORTRAN object files to be linked
- The address.F file compiled and linked for COMMON block access
- StarCallf77.h header from the build system

## Compilation

The C++ facade requires:
- C++11 or later compiler
- FORTRAN compiler for address.F
- StarCallf77.h header (provided by build system)

Example compilation (conceptual):
```bash
# Compile FORTRAN address functions
gfortran -c address.F -o address.o

# Compile C++ facade
g++ -std=c++11 -c urqmd_facade.cxx -o urqmd_facade.o

# Link with UrQMD FORTRAN library (build system specific)
```

## Verification

The implementation has been verified to:
- Compile successfully with g++ (C++11)
- Maintain the same logical structure as the FORTRAN original
- Preserve all event loop logic, collision handling, and output calls
- Properly access FORTRAN COMMON blocks through address functions

## Notes

- The facade preserves the procedural nature of the original FORTRAN code
- No object-oriented wrappers are used - this is intentional for minimal overhead
- Memory layout of COMMON blocks must match exactly between C++ and FORTRAN
- The build system provides StarCallf77.h with platform-specific F77_NAME macro
