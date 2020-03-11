test_that("R and C versions of CCMPP return the same values", {
    data(Thailand_demog)
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
    expect_equal(round(proj_r$female, 6), round(proj_c$female, 6))
    expect_equal(round(proj_r$male, 6), round(proj_c$male, 6))
})
