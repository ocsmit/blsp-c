#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "blsp_r.h"


static const R_CallMethodDef CallEntries[] = {
  {"run_blsp", (DL_FUNC) &run_blsp, 7},
  {NULL, NULL, 0}
};

void R_init_blspR(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
  R_forceSymbols(dll, TRUE);
}
