/* FILE: cohort_comp_pop_proj.c
 *
 * AUTHOR: Mark Wheldon
 *
 * DESC: C function to perform the calculations in cohort component method of
 * population projection (CCMPP). Called by the R function
 * 'ccmp_c'.
 *
 * ----------------------------------------------------------------------------
 *
 * # SYNOPSIS
 *
 * Workhorse function to compute the projected population.
 *
 * ----------------------------------------------------------------------------
 * */


/*
 * # SUB-FUNCTIONS
 *
 * */

/*
 * ## Index converter.
 *
 * */

int loc(int nr, int r, int c) {
  /* Gives location (0-based) of an element in the vector version of the element in
   * a matrix with 'nr' rows, referenced by the (1-based) array indices 'r' and
   * 'c'. */
  return (c - 1) * nr + r - 1;
}


/*
 * ## Projection Function
 *
 * */

/* Projection for a single age-group in a projection step. Adds half net
 * number of migrants, survives across projection interval, adds remaining
 * migrants. */
double projection(double out_pop_prev, double out_pop_prev_pl_1
		  ,double surv_prev
		  ,double mig_prop_prev, double mig_prop_prev_pl_1
		  ) {
  /* 'younger' and 'curr' refer to age groups. All quantities are for the
   * youngerious time period. */
  return out_pop_prev * (1 + mig_prop_prev) * surv_prev +
    mig_prop_prev_pl_1 * out_pop_prev_pl_1;
}


/*
 * # MAIN FUNCTION
 *
 * */

void ccmpp(double *out_pop_fem, double *out_pop_male
                              ,double *srb, double *fert
                              ,double *surv_fem, double *surv_male
                              ,double *surv_fem_0, double *surv_male_0
                              ,double *mig_fem, double *mig_male
                              ,double *bline_fem, double *bline_male
                              ,int *n_proj_steps /* 1-based */
                              ,int *n_age_pop
                              ,int *step_wid
    ) {

    /* INPUTS:
     *
     * - population counts (age- time- sex-specific)
     * - srb (time-specific)
     * - fertility rates (age- time-specific)
     * - survival proportions (age- time- sex-specific)
     * - survival proportions for birth to first age group (age- time- sex-specific)
     * - migration proportions (age- time- sex-specific)
     * - baseline population (age- sex-specific)
     * - number of projection steps
     * - width of the projection step and age groups
     *
     *
     * EFFECT:
     *
     * - Fills population count input vector so that it contains the projected
         population counts (age- time-specific).
     *
     * */

    /*
     * # DECLARE VARIABLES
     *
     * */

    int i, i_start, i_start_prev1, i_prev, i_prev_pl_1, y, y_prev2;
    int n_elements = *n_age_pop * *n_proj_steps;
    double fert_whole_int[n_elements];
    double mig_fem_prop_whole_int_halved[n_elements], mig_male_prop_whole_int_halved[n_elements]
        ,net_mig_fem0, net_mig_male0;
    double fert_prime, fert_next, births, srb_denom;


    /*
     * # SET-UP
     *
     * */

    /* fert and mig_prop_SEX are the average annual fertility rates and
     * migration proportions. Need to multiply by the width of the projection
     * interval to get the right values for projection. */
    for(i = 0; i < n_elements; i++) {
        fert_whole_int[i] = *step_wid * fert[i];
        mig_fem_prop_whole_int_halved[i] = *step_wid * mig_fem[i] / 2;
        mig_male_prop_whole_int_halved[i] = *step_wid * mig_male[i] / 2;
    }


    /*
     * # STEP-BY-STEP PROJECTION
     *
     * */

    for(y = 2; y <= *n_proj_steps + 1; y++) {

        /* y is the index of the year for which projection is to be computed, where
         * y = 0 is the baseline year (first "column" of 'out.pop.SEX'). */

        /* 'out_pop_SEX' is a vector converted from an R-side matrix containing
         * the population counts at baseline (first column) and projected counts
         * (subsequent columns). Projection proceeds in a sequence of
         * '*n_proj_steps' steps. Projected counts at step y are derived by
         * multiplying the counts derived in step y-1 by a function of fertility
         * rates, survival proportions and migration proportions.
         *
         * [!!!] Henceforth, I shall refer to the first column of R-side
         * matrices as 'column 1'.
         *
         * On the R-side, the population count matrix has the baseline count in
         * column 1. Projected counts are to end up in columns 2, ...,
         * *n_proj_steps. Thus, projected counts in column 1 < y <= '*n_proj_steps'
         * are a function of the counts in column y-1 and the entries in column
         * y-1 of the fert, surv and mig matrices.
         *
         * Matrix (row, column) indices need to be converted into vector indices
         * for use with the C-side versions of the R-side population count, fert
         * rate, surv prop and mig prop matrices because this is a conversion of
         * an R-program. The index 'i_start', computed using 'loc(y)', where y
         * is the index of the current projection step, indexes the first
         * element of the block of elements in the C-side vector 'out_pop_SEX'
         * into which the projected counts for this step will be placed. E.g.,
         * for step 1, 'i_start' = '*n_age_pop', which corresponds to the first
         * element of column 1 in the R-side population count matrix. In
         * general, projection step y, y = 1, ..., '*n_proj_steps', will compute
         * counts to be placed in column y + 1 of the R-side population count
         * matrix. */

        /* Find the first element of the 'column' of out_pop_SEX to be filled
         * for this value of y. Note that this column will be filled with values
         * computed from the previous column (and the corresponding columns of
         * fert, surv, mig), so we need the index of the first element of this
         * previous column as well. This is 'i_prev'. */
        i_start = loc(*n_age_pop, 1, y);

        /* Will need these a few times */
        y_prev2 = y - 2;  /* The '2' is to remind that 2 is
                           * subtracted. 'y_prev2' is used to reference
                           * quantities in the previous step but 2 is subtracted
                           * instead of 1 to make it work for 0-based
                           * indexing. */
        i_start_prev1 = i_start - *n_age_pop;

        /* NOW for the projection:
         *
         * */

        /* First age group */

        /* Compute all (female and male) births. */
        births = 0;
        /* Loop over elements of what is column y-1 of the R-side fertility rate
         * and survival proportion matrices. 'fert_prime' is the average
         * fertility rate over the projection step. It depends on the
         * fertility rates in the current and subsequent age groups, as well as
         * survival in the current age group. */
        for(i_prev = i_start_prev1; i_prev < i_start; i_prev++) {
            if(i_prev == i_start - 1) fert_next = 0;
            else fert_next = fert_whole_int[i_prev + 1];
            fert_prime = (fert_whole_int[i_prev] + fert_next * surv_fem[i_prev]) / 2;
            births = births +
                (fert_prime *
                 out_pop_fem[i_prev] * (1 + mig_fem_prop_whole_int_halved[i_prev])
                );
        }

        /* Net number of migrants in first age group are the half added at the
         * end of the projection step and not survived (this comes from the
         * matrix formulation of ccmpp---see documentation). */
        net_mig_fem0 = out_pop_fem[i_start_prev1] *
            mig_fem_prop_whole_int_halved[i_start_prev1];
        net_mig_male0 = out_pop_male[i_start_prev1] *
            mig_male_prop_whole_int_halved[i_start_prev1];

        /* Project previous period's first age group, females and males. No
         * migration applied to the births but do apply migration to existing
         * population in first age group. */
        srb_denom = (1 + srb[y_prev2]);
        out_pop_fem[i_start] = births / srb_denom * surv_fem_0[y_prev2] +
            net_mig_fem0;
        out_pop_male[i_start] = births * srb[y_prev2] / srb_denom *
            surv_male_0[y_prev2] + net_mig_male0;

        /* Project all but the first and last age groups of previous period,
         * females and males. */
        for(i = i_start + 1; i < i_start + *n_age_pop; i++) {
            i_prev = i - *n_age_pop - 1;
            i_prev_pl_1 = i_prev + 1;
            out_pop_fem[i] =
                projection(out_pop_fem[i_prev], out_pop_fem[i_prev_pl_1]
                           ,surv_fem[i_prev]
                           ,mig_fem_prop_whole_int_halved[i_prev]
                           ,mig_fem_prop_whole_int_halved[i_prev_pl_1]);
            out_pop_male[i] =
                projection(out_pop_male[i_prev], out_pop_male[i_prev_pl_1]
                           ,surv_male[i_prev]
                           ,mig_male_prop_whole_int_halved[i_prev]
                           ,mig_male_prop_whole_int_halved[i_prev_pl_1]);
        }

        /* Project previous period's last (open-ended) age group. */

        /* The projected count in the last age group, call it 'A', is
         *
         * n_{A,t+1} =
         *   n_{A-5,t}(1 + g_{A-5,t}/2)s_{A,t} + g_{A,t}/2 +
         *   n_{A,t}(1 + g_{A,t}/2)s_{A+5,t}
         *
         * The for loop computed and filled the first line. Now compute the last line. */

        /* 'i' has been left at the value i_start + *n_age_pop so we need to
         * wind it back one. 'i_prev' and i_prev_pl_1 will _not_ need to be
         * re-calculated. */
        --i;
        out_pop_fem[i] = out_pop_fem[i] +
            surv_fem[i_prev_pl_1] *
            out_pop_fem[i_prev_pl_1] * (1 + mig_fem_prop_whole_int_halved[i_prev_pl_1]);
        out_pop_male[i] = out_pop_male[i] +
            surv_male[i_prev_pl_1] *
            out_pop_male[i_prev_pl_1] * (1 + mig_male_prop_whole_int_halved[i_prev_pl_1]);
    }

    return;
}



/* ----------------------------------------------------------------------------
 * REGISTRATION
 */

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

/* .C calls */
extern void ccmpp(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *);

static const R_CMethodDef CEntries[] = {
    {"ccmpp", (DL_FUNC) &ccmpp, 15},
    {NULL, NULL, 0}
};

attribute_visible  // optional
void R_init_ccmpp(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);    
    R_forceSymbols(dll, TRUE);
}
