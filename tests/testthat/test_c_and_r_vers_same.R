test_that("R and C versions of CCMPP return the same values", {
    data(Thailand_demog)
    ## Add some fake migration to make the test worth it
    dims <- dim(Thailand_demog$thai_mig_props$female)
    Thailand_demog$thai_mig_props$female <-
        Thailand_demog$thai_mig_props$male <-
            matrix(c(0.01, 0.02, 0.03, 0.04, 0.05), nrow = dims[1], ncol = dims[2])
    proj_r <- with(Thailand_demog,
                   ccmpp_r(thai_base_pop_counts,
                      surv_props = thai_surv_props,
                      fert_rates=thai_fert_rates,
                      srb = thai_srb, mig_props = thai_mig_props))
    proj_c <- with(Thailand_demog,
                   ccmpp(thai_base_pop_counts,
                      surv_props = thai_surv_props,
                      fert_rates=thai_fert_rates,
                      srb = thai_srb, mig_props = thai_mig_props))
    proj_c_simple <- with(Thailand_demog, {
        ccmpp_c(out_pop_fem = matrix(rep(thai_base_pop_counts$female, 9), ncol = 9),
                out_pop_male = matrix(rep(thai_base_pop_counts$male, 9), ncol = 9),
                srb = thai_srb,
                fert = thai_fert_rates,
                surv_fem = thai_surv_props$female[-1,],
                surv_male = thai_surv_props$male[-1,],
                surv_fem_0 = thai_surv_props$female[1,],
                surv_male_0 = thai_surv_props$male[1,],
                mig_fem = thai_mig_props$female,
                mig_male = thai_mig_props$male,
                bline_fem = thai_base_pop_counts$female,
                bline_male = thai_base_pop_counts$male,
                n_proj_steps = ncol(thai_fert_rates), #8
                n_age_pop = nrow(thai_base_pop_counts$female), #17
                step_wid = 5
                )})
    expect_equal(round(proj_r$female, 6), round(proj_c$female, 6))
    expect_equal(round(matrix(proj_c_simple[[1]], ncol = 9), 6), round(proj_c$female, 6))
    expect_equal(round(proj_r$male, 6), round(proj_c$male, 6))
    expect_equal(round(matrix(proj_c_simple[[2]], ncol = 9), 6), round(proj_c$male, 6))
})
