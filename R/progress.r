simple_bar <- function(verbose = TRUE) {
    if (!verbose) return(NULL)
    width_display <- round(options()$width)
    bar <- paste0("#",
        paste0(rep("-", width_display), collapse = ""),
        "#")
    cat(bar, "\n")

}


message_switch <- function(verbose, type,...){
    args <- list(...)
    t <- format(Sys.time(), format = "%Y-%M-%D %H:%M:%S")
    if (verbose) {
        switch(EXPR = type,
        "kuresi" = cat(paste(t, "===> ", args$seed, "VS", args$query, "<===\n")))
    } else {
        return(NULL)
    }
}