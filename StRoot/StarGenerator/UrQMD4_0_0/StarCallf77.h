#ifndef __StarCallf77_h__
#define __StarCallf77_h__

// FORTRAN name mangling for different platforms
#if defined(__linux__) || defined(__unix__) || defined(__APPLE__)
  #define F77_NAME(name,NAME) name##_
  #define type_of_call
#else
  #define F77_NAME(name,NAME) name##_
  #define type_of_call
#endif

#endif
