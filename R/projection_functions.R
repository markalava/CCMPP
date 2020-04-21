##' Cohort-component population projection  (\R implementation)
##'
##' A pure \R implementation of the two-sex, female-dominant cohort-component method of
##' population projection (e.g., \cite{Preston et. al, 2001, Ch. 6}). This is a deterministic method for
##' projecting age-stratified population counts forward in time. It is compatible with the
##'
##' @param base_pop_counts Population counts at baseline. List with up to two components: \dQuote{female} and \dQuote{male} each a column vector of age specific counts.
##' @param surv_props List with up to two components: \dQuote{female} and \dQuote{male} each a matrix of survivorship proportions: the probability of reaching the age at the start of the interval for each projection interval. Years as columns, ages as rows. The first row should be nL0/(age_int*l0). The last row is survival for age_int years in the open interval.
##' @param fert_rates Matrix of (annual) age-specific fertility rates (NOT yet multiplied by \code{age_int}) for each projection interval. Years as columns, ages as rows.
##' @param srb Matrix of sex ratios at birth. Should be of dim 1 * \code{proj_steps}.
##' @param mig_props List with up to two components: \dQuote{female} and \dQuote{male} each a matrix of age-specific proportions for each projection interval, beginning with baseline. Years as columns, ages as rows.
##' @param proj_steps The number of time periods to project forward.
##' @param age_int Size of projection intervals (years).
##' @param label_dims Should output have dimnames set? (cosmetic).
##' @param base_year Label of baseline year (cosmetic).
##' @param first_age Label of first age group (cosmetic).
##' @param do_checks Logical; check inputs for validity?
##' @param verbose Logical; issue messages (cosmetic)?
##' @param return_list Logical. If female only projection, should a list be returned with one component called \dQuote{female}?
##' @return If \code{isTRUE(return_list)}, a list with up to two components, \dQuote{female} and \dQuote{male}, each a matrix of age-specific population counts, years as columns, ages as rows. Otherwise a matrix of age-specific counts for females only.
##' @author Mark Wheldon
##' @family CCMPP backend functions
##' @references
##' Preston, S. H., Heuveline, P., and Guillot, M. (2001), \emph{Demography: Measuring and Modeling Population Processes}, Malden, Massachusetts: Blackwell.
##' @examples
##' data("Thailand_demog")
##' with(Thailand_demog, ccmpp_r(thai_base_pop_counts,
##'                   surv_props = thai_surv_props,
##'                   fert_rates=thai_fert_rates,
##'                   srb = thai_srb, mig_props = thai_mig_props))
##' @export
ccmpp_r <- function(base_pop_counts, surv_props, fert_rates,
                    srb = matrix(1.05, ncol = proj_steps), mig_props,
                    proj_steps = ncol(fert_rates), age_int = 5,
                    label_dims = FALSE,
                    base_year = colnames(base_pop_counts[[1]]), first_age = 0,
                    do_checks = FALSE, verbose = FALSE,
                    return_list = FALSE) {

    ## -------* Housekeeping

    ## Two-sexes or one?
    is2sex <- identical(length(base_pop_counts), 2L)


    ## -------** Checks

    if(do_checks) {

        ## -------*** Check class of inputs

        ## Make sure base_pop_counts, surv_props and mig_props are lists w correct components
        if(!is2sex) {
            if(!is.list(base_pop_counts)) base_pop_counts <- list(female = as.matrix(base_pop_counts, ncol = 1))
            else if(is.null(dim(base_pop_counts[["female"]])))
                base_pop_counts <- list(female = as.matrix(base_pop_counts[["female"]], ncol = 1))
            if(!is.list(surv_props)) surv_props <- list(female = as.matrix(surv_props))
            if(!is.list(mig_props)) mig_props <- list(female = as.matrix(mig_props))
        } else {
            ## Make sure base_pop_counts is a matrix
            if(any(sapply(base_pop_counts, function(z) is.null(dim(z)))))
                base_pop_counts <- lapply(base_pop_counts, function(z) as.matrix(z, ncol = 1))
        }

        if(is2sex) {
            ## Female and male components of input lists must be matrices with
            ## matching dimensions
            invisible(lapply(list(base_pop_counts, surv_props, mig_props), FUN = function(z) {
                if(!identical(dim(z[["female"]]), dim(z[["male"]]))) stop("Male and female components must be of equal dimension within 'base_pop_counts', 'surv_props' and 'mig_props'")
            }
                             )
                      )
        }

        ## Make sure srb is a matrix
        if(identical(length(srb), 1L)) {
            srb <- matrix(srb, ncol = proj_steps)
        } else {
            if(!identical(ncol(srb), proj_steps))
                stop("!identical(ncol(srb), proj_steps)")
        }


        ## -------*** Check dimensions

        ## Number of years and ages in inputs must match
        input_dims <-
            lapply(list(fert_rates = fert_rates, surv_props = surv_props[["female"]], mig_props = mig_props[["female"]]
                        ,base_pop_counts = base_pop_counts[["female"]], srb = srb), FUN = "dim"
                   )

        input_ages <- sapply(input_dims[names(input_dims) != "srb"], "[[", 1)
        input_ages[["surv_props"]] <- input_ages[["surv_props"]] - 1
        if(!isTRUE(all.equal(input_ages - input_ages[[1]]
                             ,rep.int(0, length(input_ages))
                      ,check.attributes = FALSE
                             )))
            stop("Number of age groups in inputs not correct. Fert_Rates, mig_props, base_pop_counts must be equal, surv_props must have one more.")

        input_years <-
            sapply(input_dims[names(input_dims) != "base_pop_counts"], "[[", 2)
        if(!isTRUE(all.equal(input_years - input_years[[1]]
                             ,rep.int(0, length(input_years))
                      ,check.attributes = FALSE
                             )))
            stop("Number of years in inputs not equal")

    }


    ## Constants
    n_age_grps <- nrow(base_pop_counts[["female"]])
    n_surv_props <- nrow(surv_props[["female"]])

    ## Derive proj_steps from ncol(fert_rates)
    if(missing(proj_steps)) proj_steps <- ncol(fert_rates)

    ## Halve the migrants
    half_mig_props <- lapply(mig_props, "*", 0.5)


    ## -------* Loop over number of required time periods

    ## Initialize base_pop_counts_matrix and Leslie matrices
    if(is2sex) {
        base_pop_counts_lis <-
            list(female = matrix(0, nrow = n_age_grps, ncol = 1 + proj_steps)
                 ,male = matrix(0, nrow = n_age_grps, ncol = 1 + proj_steps))
        base_pop_counts_lis[["female"]][,1] <- base_pop_counts[["female"]]
        base_pop_counts_lis[["male"]][,1] <- base_pop_counts[["male"]]
    lesM_female <- matrix(0, nrow = n_age_grps, ncol = n_age_grps)
    lesM_male <- matrix(0, nrow = n_age_grps, ncol = n_age_grps)
    } else {
        base_pop_counts_lis <-
            list(female = matrix(0, nrow = n_age_grps, ncol = 1 + proj_steps))
        base_pop_counts_lis[["female"]][,1] <- base_pop_counts[["female"]]
    lesM_female <- matrix(0, nrow = n_age_grps, ncol = n_age_grps)
    }


    ## -------* Project in loop

    for(i in 1:proj_steps) {

        ## -------** Net migrants

        half_net_numb_mig_props <-
            mapply(FUN = function(a, b) a[,i] * b[,i] * age_int
                   ,half_mig_props, base_pop_counts_lis, SIMPLIFY = FALSE
                   )


        ## -------** Females

        ## Leslie Matrix for females: diagonal, then bottom right entry
        diag(lesM_female[-1,]) <- surv_props[["female"]][2:(n_age_grps),i]
        lesM_female[n_age_grps,n_age_grps] <- surv_props[["female"]][n_surv_props,i]

        ## First row
        k <- 1/(1+srb[i]) * surv_props[["female"]][1,i] * 0.5
        dbl_fert_rates <- age_int * (fert_rates[,i] + c(fert_rates[-1,i], 0) * surv_props[["female"]][-1,i])
        lesM_female[1,] <- k * dbl_fert_rates

        ## Project
        base_pop_counts_lis[["female"]][,i+1] <-
            lesM_female %*% (base_pop_counts_lis[["female"]][,i] +
                             half_net_numb_mig_props[["female"]]
                             ) + half_net_numb_mig_props[["female"]]


        ## -------** Males

        if(is2sex) {
            ## Leslie Matrix for males: diagonal, then bottom right entry.
            diag(lesM_male[-1,]) <- surv_props[["male"]][2:(n_age_grps),i]
            lesM_male[n_age_grps,n_age_grps] <- surv_props[["male"]][n_surv_props,i]
            base_pop_counts_lis[["male"]][,i+1] <-
                lesM_male %*% (base_pop_counts_lis[["male"]][,i] + half_net_numb_mig_props[["male"]]) +
                    half_net_numb_mig_props[["male"]]

            ## Add male births. Migration of mothers in this interval already
            ## taken care of in previous step.
            male_k <- srb[,i] / surv_props[["female"]][1,i] * surv_props[["male"]][1,i]
            base_pop_counts_lis[["male"]][1,i+1] <- base_pop_counts_lis[["male"]][1,i+1] +
                male_k * (base_pop_counts_lis[["female"]][1,i+1] - half_net_numb_mig_props[["female"]][1])
        }

    }


    ## -------* OUTPUT

    ## Add dim names
    if(label_dims) {
        ages <- seq(from = first_age
                    ,to = age_int*(nrow(as.matrix(base_pop_counts_lis[["female"]]))-1)
                    ,by = age_int)
        years <- (0:proj_steps) * age_int + as.numeric(base_year)
        for(i in 1:length(base_pop_counts_lis))
            dimnames(base_pop_counts_lis[[i]]) <- list(age = ages, year = years)
    }


    ## Return
    if(!is2sex && !return_list) return(base_pop_counts_lis[[1]])
    else return(base_pop_counts_lis)

}


##' Cohort-component population projection
##'
##' Female-dominant cohort-component method of
##' population projection (e.g., \cite{Preston et. al, 2001, Ch. 6}). This is a deterministic method for
##' projecting age-stratified population counts forward in time. It calls compiled code to do the actual projection.
##'
##' @seealso \code{\link{ccmpp_c}} which is a bare wrapper for the same underlying _C_ function, and \code{\link{ccmpp_r}} which is an implementation in pure \R.
##'
##' @param n_age_grps Number of age groups
##' @inheritParams ccmpp_r
##' @return If \code{isTRUE(return_list)}, a list with up to two components, \dQuote{female} and \dQuote{male}, each a matrix of age-specific population counts, years as columns, ages as rows. Otherwise a matrix of age-specific counts for females only.
##' @author Mark Wheldon
##' @references
##' Preston, S. H., Heuveline, P., and Guillot, M. (2001), \emph{Demography: Measuring and Modeling Population Processes}, Malden, Massachusetts: Blackwell.
##' @examples
##' data("Thailand_demog")
##' with(Thailand_demog, ccmpp(thai_base_pop_counts,
##'                   surv_props = thai_surv_props,
##'                   fert_rates=thai_fert_rates,
##'                   srb = thai_srb, mig_props = thai_mig_props))
##' @useDynLib ccmpp, .registration = TRUE, .fixes = "C_"
##' @export
ccmpp <- function(base_pop_counts, surv_props, fert_rates,
                    srb = matrix(1.05, ncol = proj_steps), mig_props,
                    proj_steps = ncol(fert_rates),
                    n_age_grps = length(base_pop_counts$female), age_int = 5,
                    label_dims = FALSE,
                    base_year = colnames(base_pop_counts[[1]]), first_age = 0,
                    do_checks = FALSE, verbose = FALSE,
                    return_list = FALSE) {


    ## -------* Housekeeping

    ## -------** Checks

    if(do_checks) {

        ## -------*** Check class of inputs

        ## Make sure pop is a matrix
            if(any(sapply(base_pop_counts, function(z) is.null(dim(z)))))
                base_pop_counts <- lapply(base_pop_counts, function(z) as.matrix(z, ncol = 1))

            ## Female and male components of input lists must be matrices with
            ## matching dimensions
            invisible(lapply(list(base_pop_counts, surv_props, mig_props), FUN = function(z) {
                if(!identical(dim(z[["female"]]), dim(z[["male"]]))) stop("Male and female components must be of equal dimension within 'base_pop_counts', 'surv_props' and 'mig_props'")
            }
                             )
                      )

        ## Make sure srb is a matrix
        if(identical(length(srb), 1L)) srb <- matrix(srb, ncol = proj_steps)
        else if(!identical(ncol(srb), proj_steps)) stop("!identical(ncol(srb), proj_steps)")


        ## -------*** Check dimensions

        ## Number of years and ages in inputs must match
        input_dims <-
            lapply(list(fert_rates = fert_rates, surv_props = surv_props[["female"]], mig_props = mig_props[["female"]]
                        ,base_pop_counts = base_pop_counts[["female"]], srb = srb), FUN = "dim"
                   )

        input_ages <- sapply(input_dims[names(input_dims) != "srb"], "[[", 1)
        input_ages[["surv_props"]] <- input_ages[["surv_props"]] - 1
        if(!all.equal(input_ages - input_ages[[1]], rep.int(0, length(input_ages))
                      ,check.attributes = FALSE)) stop("Number of age groups in inputs not equal")

        input_years <-
            sapply(input_dims[names(input_dims) != "base_pop_counts"], "[[", 2)
        if(!all.equal(input_years - input_years[[1]], rep.int(0, length(input_years))
                      ,check.attributes = FALSE)) stop("Number of years in inputs not equal")

    }


    ## Constants
    if(missing(proj_steps)) proj_steps <- ncol(fert_rates)

    ## Initialize pop_matrix and Leslie matrices
    n_elements <- n_age_grps * (1 + proj_steps)
    base_pop_counts_lis_fem <- base_pop_counts_lis_male <- vector("double", length = n_elements)

    base_pop_counts_lis_fem[1:n_age_grps] <- base_pop_counts$fem
    base_pop_counts_lis_male[1:n_age_grps] <- base_pop_counts$male


    ## -------* Call C program

    cfun <- .C(C_ccmpp
               ,out_pop_fem = as.double(base_pop_counts_lis_fem)
               ,out_pop_male = as.double(base_pop_counts_lis_male)
               ,srb = as.double(srb)
               ,fert_rates = as.double(fert_rates)
               ,surv_fem = as.double(surv_props$female[-1,])
               ,surv_male = as.double(surv_props$male[-1,])
               ,surv_fem_0 = as.double(surv_props$female[1,])
               ,surv_male_0 = as.double(surv_props$male[1,])
               ,mig_fem = as.double(mig_props$female)
               ,mig_male = as.double(mig_props$male)
               ,bline_fem = as.double(base_pop_counts$female)
               ,bline_male = as.double(base_pop_counts$male)
               ,n_proj_steps = as.integer(proj_steps)
               ,n_age_pop = as.integer(n_age_grps)
               ,step_wid = as.integer(age_int)
               )

    base_pop_counts_lis <- list(female = matrix(cfun[[1]], nrow = n_age_grps)
                    ,male = matrix(cfun[[2]], nrow = n_age_grps)
                    )


    ## -------* OUTPUT

    ## Add dim names
    if(label_dims) {
        ages <- seq(from = first_age
                    ,to = age_int*(nrow(as.matrix(base_pop_counts_lis[["female"]]))-1)
                    ,by = age_int)
        years <- (0:proj_steps) * age_int + as.numeric(base_year)
        for(i in 1:length(base_pop_counts_lis))
            dimnames(base_pop_counts_lis[[i]]) <- list(age = ages, year = years)
    }

    return(base_pop_counts_lis)
}


##' \R wrapper for the _C_ function \code{C_ccmpp}.
##'
##' This is a simple wrapper for the _C_ function of the same name
##' which performs population projection by the cohort component
##' method (see the documentation for \code{\link{ccmpp}}). The _C_
##' function can be loaded for use in _C_ source code via
##' \code{#include ccmppAPI}. The arguments are documented here for
##' both functions. Argument types for _C_ are given in parentheses in
##' the argument list.
##'
##' The _C_ code is accessed via \code{\link{.C}} hence all arguments
##' are passed via pointers. In particular, the first two arguments
##' are containers for the output. They must be vectors of length
##' \code{n_age_pop} which can be coerced via
##' \code{\link{as.double}}. The values of the elements are
##' unimportant as they will be overwritten with the output; e.g., you
##' could supply the same value to \code{out_pop_fem} as
##' \code{bline_fem} (and similarly for the male versions), as in the
##' example.
##'
##' @seealso \code{\link{ccmpp}}, a more friendly wrapper for
##'     the same underlying _C_ function that takes arguments as \R
##'     lists and does some optional input checks.
##'
##' @param out_pop_fem (\code{double *out_pop_fem}) In \R this is an
##'     \code{n_age_pop} * (1 + \code{proj_steps}) matrix with age
##'     groups as rows and time periods as columns. In _C_ it is an
##'     array of length \code{n_age_pop} * (1 +
##'     \code{proj_steps}). Its first column (in \R), or first
##'     \code{n_age_pop} entries (in _C_) contain \code{bline_fem},
##'     the age-specific female population counts in the year that
##'     marks the start of the projection (the \dQuote{baseline}). The
##'     rest of the entries are arbitrary; they serve as a repository
##'     for the results which will overrite whatever entries are
##'     supplied. See \dQuote{Details}.
##' @param out_pop_male (\code{double *out_pop_male}) Same as
##'     \code{out_pop_fem} but for males.
##' @param srb (\code{double *srb}) Vector of length
##'     \code{n_proj_steps} containing the average sex ratios at
##'     birth.
##' @param fert (\code{double *fert}) Average annual age-specific
##'     fertility rates for each of the \code{n_proj_steps}. In \R
##'     this can be supplied as a matrix with age groups as rows and
##'     time periods as columns. In _C_ it is an array of length
##'     \code{n_age_pop} * \code{proj_steps} formed by concatenating
##'     the age-specific fertility rates for each projection period.
##' @param surv_fem (\code{double *surv_fem}) Age-specific survival
##'     proportions for females for ages 5 to the open ended age
##'     group. Like \code{fert}, in \R it can be supplied as an
##'     \code{n_age_pop} by \code{n_proj_steps} matrix with age groups
##'     as rows. If using the _C_ version it should be an array with
##'     \code{n_age_pop} * \code{n_proj_steps} elements.
##' @param surv_male (\code{double *surv_male}) Same as
##'     \code{surv_fem} but for males.
##' @param surv_fem_0 (\code{double *surv_fem_0}) Vector of lenght
##'     \code{n_proj_steps} containing survival proportions of females
##'     for the period birth to the end of the first projection
##'     interval, for each period in the projection.
##' @param surv_male_0 (\code{double *surv_male_0}) Same as
##'     \code{surv_fem_0} but for males.
##' @param mig_fem (\code{double *mig_fem}) A matrix (in \R) or vector
##'     (in _C_) with the same dimensions as for \code{fert} but
##'     containing average annual net migrations for females as a
##'     proportion of the receiving population, by age.
##' @param mig_male (\code{double *mig_male}) Same as \code{mig_fem}
##'     but for males.
##' @param bline_fem (\code{double *bline_fem}) A vector of
##'     age-specific female population counts in the year that marks
##'     the start of the projection (the \dQuote{baseline}).
##' @param bline_male (\code{double *bline_male}) Same as
##'     \code{bline_fem} but for males.
##' @param n_proj_steps (\code{int *n_proj_steps}) The number of steps
##'     to project forward. This must conform to the dimensions of the
##'     other inputs; see where this argument is mentioned in their
##'     descriptions. \emph{No checking is done} to make sure this
##'     condition is satisfied.
##' @param n_age_pop (\code{int *n_age_pop}) The number of age groups
##'     in the projection. As for \code{n_proj_steps}, this must
##'     conform to the dimensions of the other inputs; see where this
##'     argument is mentioned in their descriptions. \emph{No checking
##'     is done} to make sure this condition is satisfied.
##' @param step_wid (\code{int *step_wid}) The scale of the projection
##'     intervals, e.g., in years.
##'
##' @examples
##' data("Thailand_demog")
##' raw <-
##'   with(Thailand_demog, {
##'   ccmpp_c(out_pop_fem = matrix(rep(thai_base_pop_counts$female, 9), ncol = 9),
##'           out_pop_male = matrix(rep(thai_base_pop_counts$male, 9), ncol = 9),
##'           srb = thai_srb,
##'           fert = thai_fert_rates,
##'           surv_fem = thai_surv_props$female[-1,],
##'           surv_male = thai_surv_props$male[-1,],
##'           surv_fem_0 = thai_surv_props$female[1,],
##'           surv_male_0 = thai_surv_props$male[1,],
##'           mig_fem = thai_mig_props$female,
##'           mig_male = thai_mig_props$male,
##'           bline_fem = thai_base_pop_counts$female,
##'           bline_male = thai_base_pop_counts$male,
##'           n_proj_steps = ncol(thai_fert_rates), #8
##'           n_age_pop = nrow(thai_base_pop_counts$female), #17
##'           step_wid = 5
##'   )})
##'
##' (nicer <- list(female = matrix(raw[[1]], ncol = 9),
##'               male = matrix(raw[[2]], ncol = 9)))
##'
##' @return In \R, a two-element list containig the first two
##'     arguments with their (\code{n_age_pop} + 1), (\code{n_age_pop}
##'     + 2), ..., elements replaced with the projected counts for
##'     females and males, respectiely. The _C_ function returns
##'     nothing (\code{void}) but has the side effect of similarly
##'     replacing element in the first two arguments.
##' @author Mark C Wheldon
##' @family CCMPP backend functions
##' @useDynLib ccmpp, .registration = TRUE, .fixes = "C_"
##' @export
ccmpp_c <- function(out_pop_fem, out_pop_male, srb, fert, surv_fem, surv_male,
                    surv_fem_0, surv_male_0, mig_fem, mig_male,
                    bline_fem, bline_male, n_proj_steps,
                    n_age_pop, step_wid) {
    .C(C_ccmpp,
       as.double(out_pop_fem),
       as.double(out_pop_male),
       as.double(srb),
       as.double(fert),
       as.double(surv_fem),
       as.double(surv_male),
       as.double(surv_fem_0),
       as.double(surv_male_0),
       as.double(mig_fem),
       as.double(mig_male),
       as.double(bline_fem),
       as.double(bline_male),
       as.integer(n_proj_steps),
       as.integer(n_age_pop),
       as.integer(step_wid))[1:2]
}

