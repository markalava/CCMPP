#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>

void C_ccmpp(double *out_pop_fem,
	   double *out_pop_male,
	   double *srb, double *fert,
	   double *surv_fem,
	   double *surv_male,
	   double *surv_fem_0,
	   double *surv_male_0,
	   double *mig_fem,
	   double *mig_male,
	   double *bline_fem,
	   double *bline_male,
	   int *n_proj_steps, /* 1-based */
	   int *n_age_pop,
	   int *step_wid) {
  static void(*fun)(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *) = NULL;
  if (fun == NULL)
    fun = (void(*)(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *)) R_GetCCallable("ccmpp", "C_ccmpp");
  return fun(*out_pop_fem,
	    *out_pop_male,
	    *srb,  *fert,
	    *surv_fem,
	    *surv_male,
	    *surv_fem_0,
	    *surv_male_0,
	    *mig_fem,
	    *mig_male,
	    *bline_fem,
	    *bline_male,
	   *n_proj_steps, /* 1-based */
	   *n_age_pop,
	   *step_wid);
}
