#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(xfourhr)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xpwe)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xpwecxpwu)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xpwecxpwuforvar)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xpwefv2)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xpwefvplus)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xpwu)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xqpwe)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xqpwu)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xrmsth)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xspf)(void *, void *, void *, void *);
extern void F77_NAME(xwlrcal)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(xwlrutil)(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"xfourhr",         (DL_FUNC) &F77_NAME(xfourhr),         10},
  {"xpwe",            (DL_FUNC) &F77_NAME(xpwe),             6},
  {"xpwecxpwu",       (DL_FUNC) &F77_NAME(xpwecxpwu),       19},
  {"xpwecxpwuforvar", (DL_FUNC) &F77_NAME(xpwecxpwuforvar), 20},
  {"xpwefv2",         (DL_FUNC) &F77_NAME(xpwefv2),          8},
  {"xpwefvplus",      (DL_FUNC) &F77_NAME(xpwefvplus),      14},
  {"xpwu",            (DL_FUNC) &F77_NAME(xpwu),             6},
  {"xqpwe",           (DL_FUNC) &F77_NAME(xqpwe),            6},
  {"xqpwu",           (DL_FUNC) &F77_NAME(xqpwu),            6},
  {"xrmsth",          (DL_FUNC) &F77_NAME(xrmsth),          10},
  {"xspf",            (DL_FUNC) &F77_NAME(xspf),             4},
  {"xwlrcal",         (DL_FUNC) &F77_NAME(xwlrcal),         14},
  {"xwlrutil",        (DL_FUNC) &F77_NAME(xwlrutil),         8},
  {NULL, NULL, 0}
};

void R_init_PWEALL(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
