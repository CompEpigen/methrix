#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

// External function declaration
extern SEXP read_modkit_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// Define the C methods that can be called from R
static const R_CallMethodDef CallEntries[] = {
  {"read_modkit_c", (DL_FUNC) &read_modkit_c, 9},
  {NULL, NULL, 0}
};

// Initialize the dynamic library
void R_init_methrix(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  
  // Optional: Print debug message
  // Rprintf("Registered read_modkit_c function\n");
}