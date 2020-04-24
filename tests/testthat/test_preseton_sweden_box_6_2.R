test_that("All functions reproduce Preston et al. (2011) Box 6.2", {

    ## -------* INPUTS

    ## -------** C Versions

    step_wid <- 5
    srb_c <- 1.05
    fert <- c(0, 0, 0, 0.012, 0.0908, 0.1499, 0.1125, 0.0441, 0.0074, 3e-04,
              0, 0, 0, 0, 0, 0, 0)
    surv_fem <- surv_male <-
        c(0.999298474130982, 0.999523271204374, 0.999255384875458, 0.998733211018043,
          0.998519868845054, 0.998071361638878, 0.997059972117285, 0.995306031333712,
          0.992182641549444, 0.987776734283712, 0.981527675199983, 0.970465422882409,
          0.952606955959568, 0.923471809204568, 0.86986034386444, 0.774955902248557,
          0.517928194893993)
    surv_fem_0 <- surv_male_0 <- 0.994974
    bline_fem <- bline_male <-
        c(293395, 248369, 240012, 261346, 285209, 314388, 281290, 286923, 304108,
          324946, 247613, 211351, 215140, 221764, 223506, 183654, 254414)
    mig_fem <- mig_male <-
        c(0.00466265614615109, 0.00334180191569801, 0.0028040264653434, 0.00403296778982651,
          0.00647945892310551, 0.00523556878761276, 0.00388922464360624, 0.00219919630005263,
          0.00116406013653044, 0.000686267872200304, 0.000868290437093367, 0.000799617697574178,
          0.000599609556567816, 0.000477985606320232, 0.000416096212182223, 0.000326701297004149,
          0.000334101110787928)
    n_proj_steps <- 1
    n_age_pop <- length(fert)
    out_pop_fem <- out_pop_male <-
        matrix(c(bline_fem, rep(0, n_age_pop)), ncol = 2,
               dimnames = list(seq(from = 0, by = 5, length = n_age_pop), NULL))

    ## -------** R Versions

    ages_char <- seq(from = 0, length = n_age_pop, by = 5)
    base_pop_counts <-
    list(female = matrix(bline_fem, dimnames = list(ages_char, NULL)),
         male = matrix(bline_male, dimnames = list(ages_char, NULL)))
    surv_props <-
    list(female = matrix(c(surv_fem_0, surv_fem),
                         dimnames = list(c(ages_char, "85"), NULL)),
         male = matrix(c(surv_male_0, surv_male),
                                dimnames = list(c(ages_char, "85"), NULL)))
    fert_rates <- matrix(fert, dimnames = list(ages_char, NULL))
    srb_R <- matrix(1.05, ncol = n_proj_steps)
    mig_props <- list(female = matrix(mig_fem, dimnames = list(ages_char, NULL)),
                  male = matrix(mig_male, dimnames = list(ages_char, NULL)))
    proj_steps <- n_proj_steps
    age_int <- step_wid

    ## -------* Outputs to Check Against

    fem_pop_check <-
        c(302392, 298682, 252007, 244150, 268267, 293515, 320624,
          284767, 288031, 303166, 322062, 243989, 205841, 205516,
          205270, 194771, 142565 + 131966)
    male_pop_check <-
        c(317340, 298682, 252007, 244150, 268267, 293515, 320624,
          284767, 288031, 303166, 322062, 243989, 205841, 205516,
          205270, 194771, 142565 + 131966)

    ## -------* TESTs

    r_proj <- ccmpp_r(base_pop_counts = base_pop_counts,
                 surv_props = surv_props, fert_rates = fert_rates,
                 srb = srb_R, mig_props = mig_props,
                 proj_steps = proj_steps,
                 label_dims = FALSE,
                 base_year = 1993, first_age = 0,
                 do_checks = TRUE, verbose = FALSE,
                 return_list = FALSE
            )

    expect_equal(round(r_proj$female[,2]), fem_pop_check)
    expect_equal(round(r_proj$male[,2]), male_pop_check)

    c_proj <- ccmpp_c(out_pop_fem = out_pop_fem, out_pop_male = out_pop_male,
                 srb = srb_c, fert = fert,
                 surv_fem = surv_fem, surv_male = surv_male,
                 surv_fem_0 = surv_fem_0, surv_male_0 = surv_male_0,
                 mig_fem = mig_fem, mig_male = mig_male,
                 bline_fem = bline_fem, bline_male = bline_male,
                 n_proj_steps = n_proj_steps, n_age_pop = n_age_pop,
                 step_wid = step_wid)

    expect_equal(round(c_proj[[1]][(n_age_pop + 1):(2 * n_age_pop)]), fem_pop_check)
    expect_equal(round(c_proj[[2]][(n_age_pop + 1):(2 * n_age_pop)]), male_pop_check)
    })
