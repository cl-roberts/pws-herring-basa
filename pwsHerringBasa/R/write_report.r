#' Write model report 
#'
#' Utility function for writing model estimates of estimated parameters and 
#' derived quantities.
#'
#' @param obj A list object returned by `MakeADFun`
#' @param par A list of parameter values to generate report
#' @param file File path to write report
#' @param llik Character vector of names of log-likelihood components
#' @param derived Character vector of names of derived quantities
#' @param survey Character vector of names of survey quantities
#' @param forecast Character vector of names of forecast quantities
#'

write_report <- function(obj, par, file, llik = NULL, derived = NULL, 
                         survey = NULL, forecast = NULL) {

    # local vars
    nyr <- obj$env$data$nyr

    # report sections
    report <- obj$report(par)
    llik_report <- report[names(report) %in% llik]
    derived_report <- report[names(report) %in% derived]
    survey_report <- report[names(report) %in% survey]
    forecast_report <- report[names(report) %in% forecast]
    
    # local functions

    # this utility function concatenates lists of objects to write and labels 
    # them with strings for writing report file
    reporter <- function(x) {

        for (i in 1:length(x)) {

            if (i == 1) {
                out <- c(paste("#", names(x)[i]), x[i])
            } else if (i > 1) {
                out <- c(out, paste("\n#", names(x)[i]), x[i])
            }

        }

        return(out)
    
    }

    # this utility function makes a list from a vector and collates all like-named
    # vector elements

    collate_list <- function(x) {

        list_names <- unique(names(x))

        out <- vector(mode = "list")

        for (i in list_names) {
            out[[i]] <- x[names(x) %in% i]
        }

        return(out)
    
    }

    # this utility function writes different sections of the report

    write_section <- function(x, file) {
            y <- as.matrix(x)
            if (nrow(y) < nyr) y <- t(y)
            ncolumns <- ncol(y)
            write(t(y), file = file, append = TRUE, ncolumns = ncolumns)
        }

    # write report file
    write(
        "# -------- PWS_ASA TMB Report File -------- #", 
        file = file
    )

    write(
        "\n\n# General information ---- \n", 
        file = file, append = TRUE
    )

    write(
        paste("# Number of parameters:", length(par)), 
        file = file, append = TRUE
    )

    write(
        paste("# Objective function value:", obj$fn(par)), 
        file = file, append = TRUE
    )

    write(
        paste("# Maximum gradient:", max(abs(obj$gr(par)))), 
        file = file, append = TRUE
    )

    write(
        "\n\n# Estimated parameters ---- \n", 
        file = file, append = TRUE
    )

    lapply(
        reporter(collate_list(par)),
        FUN = \(x) write_section(x, file = file)
    )

    if (!is.null(llik)) {

        write(
            "\n\n# Likelihood components ---- \n", 
            file = file, append = TRUE
        )

        lapply(
            reporter(llik_report), 
            FUN = \(x) write_section(x, file = file)
        )
        
    }

    if (!is.null(derived)) {

        write(
            "\n\n# Derived quantities ---- \n", 
            file = file, append = TRUE
        )

        lapply(
            reporter(derived_report), 
            FUN = \(x) write_section(x, file = file)
        )

    }

    if (!is.null(survey)) {

        write(
            "\n\n# Survey quantities ---- \n", 
            file = file, append = TRUE
        )

        lapply(
            reporter(survey_report), 
            FUN = \(x) write_section(x, file = file)
        )

    }

    if (!is.null(forecast)) {

        write(
            "\n\n# Forecast quantities ---- \n", 
            file = file, append = TRUE
        )

        lapply(
            reporter(forecast_report), 
            FUN = \(x) write_section(x, file = file)
        )
    
    }

    return(0)

}
