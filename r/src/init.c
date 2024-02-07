#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "blsp_r.h"



static const R_CallMethodDef CallEntries[] = {
  {"generate_data", (DL_FUNC) &generate_data, 3},
  {NULL, NULL, 0}
};

void R_init_blspR(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
  R_forceSymbols(dll, TRUE);
}
