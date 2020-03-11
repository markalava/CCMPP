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
##' @param n_age_grps Number of age groups
##' @inheritParams ccmpp_r
##' @return If \code{isTRUE(return_list)}, a list with up to two components, \dQuote{female} and \dQuote{male}, each a matrix of age-specific population counts, years as columns, ages as rows. Otherwise a matrix of age-specific counts for females only.
##' @author Mark Wheldon
##' @family CCMPP backend functions
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
               ,out.pop.fem = base_pop_counts_lis_fem
               ,out.pop.male = base_pop_counts_lis_male
               ,srb = as.double(srb)
               ,fert_rates = as.double(fert_rates)
               ,surv.fem = as.double(surv_props$female[-1,])
               ,surv.male = as.double(surv_props$male[-1,])
               ,surv.fem.0 = as.double(surv_props$female[1,])
               ,surv.male.0 = as.double(surv_props$male[1,])
               ,mig.fem = as.double(mig_props$female)
               ,mig.male = as.double(mig_props$male)
               ,bline.fem = as.double(base_pop_counts$female)
               ,bline.male = as.double(base_pop_counts$male)
               ,n.proj.steps = as.integer(proj_steps)
               ,n.age.pop = as.integer(n_age_grps)
               ,step.wid = as.integer(age_int)
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
