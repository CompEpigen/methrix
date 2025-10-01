#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

// External function declarations
extern SEXP read_modkit_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP read_modkit_v2_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP add_context_c(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP read_modkit_v2_init_chunks_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP read_modkit_v2_process_chunk_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP read_modkit_v2_cleanup_c(SEXP, SEXP);
extern SEXP read_modkit_v2_process_interval_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP read_modkit_v2_create_intervals_c(SEXP, SEXP, SEXP, SEXP);

// Define the C methods that can be called from R
static const R_CallMethodDef CallEntries[] = {
  {"read_modkit_c", (DL_FUNC) &read_modkit_c, 10},
  {"read_modkit_v2_c", (DL_FUNC) &read_modkit_v2_c, 9},
  {"add_context_c", (DL_FUNC) &add_context_c, 5},
  {"read_modkit_v2_init_chunks_c", (DL_FUNC) &read_modkit_v2_init_chunks_c, 9},
  {"read_modkit_v2_process_chunk_c", (DL_FUNC) &read_modkit_v2_process_chunk_c, 9},
  {"read_modkit_v2_cleanup_c", (DL_FUNC) &read_modkit_v2_cleanup_c, 2},
  {"read_modkit_v2_process_interval_c", (DL_FUNC) &read_modkit_v2_process_interval_c, 9},
  {"read_modkit_v2_create_intervals_c", (DL_FUNC) &read_modkit_v2_create_intervals_c, 4},
  {NULL, NULL, 0}
};

// Initialize the dynamic library
void R_init_methrix(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  
  // Optional: Print debug message
  // Rprintf("Registered read_modkit_c function\n");
}