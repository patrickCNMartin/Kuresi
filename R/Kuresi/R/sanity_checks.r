
check_comp_genes <- function(signal, vesalius_assay, norm_method, method) {
    method <- vesalius:::check_deg_method(method)
    counts <- vesalius:::check_norm(vesalius_assay, norm_method, method, FALSE)
    signal <- vesalius:::check_genes(signal, counts)
    return(signal)
}