#' Retrieve Implemented Test Problems for the SOLNP Suite
#'
#' Returns a list (or a single object) of implemented test problems corresponding to a selected suite.
#' Problem functions must follow the naming convention \sQuote{<problem_name>_problem} and return a list
#' describing the optimization problem (e.g., objective, constraints, bounds).
#'
#' @param suite Character. The test suite to draw from. Must be one of \dQuote{Hock-Schittkowski} or \dQuote{Other}.
#'              Default is \dQuote{Hock-Schittkowski}.
#' @param number Integer or vector of integers. One or more problem numbers to retrieve. Ignored if return_all = TRUE.
#' @param return_all Logical. If TRUE, returns all implemented problems in the specified suite. Default is FALSE.
#'
#' @return If one problem is requested and implemented, the evaluated problem object is returned directly.
#'         Otherwise, an unnamed list of evaluated problem objects is returned.
#'
#' @details
#' - Problems are matched by number within the selected suite, using the table from [solnp_problems_table()].
#' - If a requested problem is valid but not yet implemented (i.e., the corresponding function does not exist),
#'   a message will inform the user.
#' - If a problem number exceeds the allowable range (e.g., > 306 for Hock-Schittkowski), an error is raised.
#'
#' @examples
#' \dontrun{
#' # Retrieve a single HS problem
#' prob <- solnp_problem_suite(number = 1)
#'
#' # Retrieve multiple HS problems
#' probs <- solnp_problem_suite(number = c(1, 2, 3))
#'
#' # Retrieve all implemented problems in "Other" suite
#' other_prob <- solnp_problem_suite(suite = "Other", number = 1)
#' }
#'
#' @seealso [solnp_problems_table()]
#' @export
solnp_problem_suite <- function(suite = "Hock-Schittkowski", number = 1, return_all = FALSE)
{
    suite <- match.arg(suite[1], c("Hock-Schittkowski","Other"))
    v_problems <- solnp_problems_table()
    v_problems <- v_problems[v_problems$Suite == suite,]
    if (!return_all) {
        number <- unique(as.integer(number))
        if (any(number < 1)) stop("Problem number must be strictly positive.")
        if (suite == "Hock-Schittkowski" & any(number > 306)) stop("Some requested problem numbers do not exist (maximum is 306 for Hock-Schittkowski).")
        known <- v_problems$Number
        not_implemented <- setdiff(number, known)
        implemented <- intersect(number, known)

        if (length(not_implemented) > 0) {
            message("The following problems are not yet implemented: ",
                    paste(not_implemented, collapse = ", "))
        }
        funs <- lapply(implemented, function(n) {
            fname <- paste0(v_problems$Problem[v_problems$Number == n], "_problem")
            if (!exists(fname, mode = "function", inherits = TRUE)) {
                warning(sprintf("Function '%s' not found in the environment.", fname))
                return(NULL)
            }
            get(fname, mode = "function", inherits = TRUE)()
        })
        funs <- funs[!sapply(funs, is.null)]
        if (length(funs) == 1) {
            return(funs[[1]])
        } else {
            return(unname(funs))
        }
    } else {
        # return_all = TRUE: Return all implemented hsXX_problem functions
        all_funs <- lapply(seq_len(nrow(v_problems)), function(i) {
            fname <- paste0(v_problems$Problem[i], "_problem")
            if (exists(fname, mode = "function", inherits = TRUE)) {
                get(fname, mode = "function", inherits = TRUE)()
            } else {
                NULL
            }
        })
        all_funs <- all_funs[!sapply(all_funs, is.null)]
        return(unname(all_funs))
    }
}

#' List of Valid Test Problems for the SOLNP Suite
#'
#' Returns a data.frame of known and registered test problems used with the SOLNP solver.
#' The list includes problems from the Hock-Schittkowski suite as well as a selection
#' of other classic optimization problems.
#'
#' @return A data.frame with the following columns:
#' \describe{
#'   \item{Suite}{A character string indicating the suite the problem belongs to.
#'                One of \dQuote{Hock-Schittkowski} or \dQuote{Other}.}
#'   \item{Problem}{The base name of the problem function (without the \sQuote{_problem} suffix).}
#'   \item{Number}{An integer identifier used to index or request problems programmatically.}
#' }
#'
#' @details
#' - All problem functions are expected to follow the naming convention
#'   \sQuote{<Problem>_problem} (e.g., \sQuote{hs01_problem}).
#' - For Hock-Schittkowski problems, numbers range from 1 to 50, with a few
#'   selected extras (e.g., 110, 118, 119).
#' - The \dQuote{Other} suite includes named problems like \sQuote{box},
#'   \sQuote{alkylation}, \sQuote{entropy}, \sQuote{garch}, etc.,
#' and are numbered sequentially.
#'
#' @examples
#' # View all known problems
#' solnp_problems_table()
#'
#' # Filter only HS problems
#' subset(solnp_problems_table(), Suite == "Hock-Schittkowski")
#'
#' @seealso [solnp_problem_suite()]
#' @export
solnp_problems_table <- function()
{
    hs_problems <- c("hs01","hs02","hs03","hs04","hs05","hs06","hs07",
                     "hs08","hs09", paste0("hs",10:56),"hs110","hs118","hs119")
    hs_number <- c(1:56,110,118,119)
    other_problems <- sort(c("alkylation","wright4","wright9","box","entropy","rosen_suzuki","powell",
                             "himmelblau5","garch"))
    other_number <- c(1:length(other_problems))
    rbind(data.frame("Suite" = "Hock-Schittkowski", "Problem" = hs_problems, "Number" = hs_number),
          data.frame("Suite" = "Other", "Problem" = other_problems, "Number" = other_number))
}
