#' Split SingleCellExperiment by colData variable
#'
#' @param x SingleCellExperiment object
#' @param f colData variable as a string
#'
#' @return a list of singlecellexperiments name by colData value
#' @export
#'
#' @examples
#' 
#' 
#' splitByCol(small_example_dataset, "batch")
splitByCol <- function(x, f = "batch") {
    f <- colData(x)[[f]]

    i <- split(seq_along(f), f)

    v <- vector(mode = "list", length = length(i))

    names(v) <- names(i)

    for (n in names(i)) {
        v[[n]] <- x[, i[[n]]]
    }

    return(v)
}


#' Merge Small SingleCellExperiment Objects
#'
#' @param ... two or more singlecell objects
#' @param k.filter minimum cell number for integration
#'
#' @return a SingleCellExperiment object
merge_small_objects <- function(..., k.filter = 50) {
    object_list <- list(...)

    # check if any singlecell objects are too small and if so merge with the first singlecell objects
    object_dims <- map(object_list, dim) %>%
        map_lgl(~ .x[[2]] < k.filter)

    small_objects <- object_list[object_dims]

    object_list <- object_list[!object_dims]

    object_list[[1]] <- reduce(c(small_objects, object_list[[1]]), correctExperiments, PARAM = NoCorrectParam())

    return(object_list)
}


#' Batch Correct Multiple Single Cell Objects
#'
#' @param object_list List of two or more SingleCellExperiment objects
#' @param organism human or mouse
#' @param ... extra args passed to object_reduce_dimensions
#'
#' @return an integrated SingleCellExperiment object
object_integrate <- function(object_list, organism = "human", ...) {
    # drop 'colData' fields with same name as 'batchCorrect' output
    object_list <- map(object_list, ~ {
        colData(.x)[["batch"]] <- NULL
        return(.x)
    })

    geneCorrected <- correctExperiments(object_list)
    mainExpName(geneCorrected) <- "integrated"

    geneMerged <- correctExperiments(object_list)
    altExp(geneCorrected, "gene") <- geneMerged

    alt_exp_names <- map(object_list, altExpNames)

    if (all(map_lgl(alt_exp_names, ~ {
        length(.x) > 0 & "transcript" %in% .x
    }))) {
        # drop 'colData' fields with same name as 'batchCorrect' output
        object_list <- map(object_list, ~ {
            colData(.x)[["batch"]] <- NULL
            return(.x)
        })

        transcriptBatches <- map(object_list, swapAltExp, "transcript")
        transcriptMerged <- correctExperiments(transcriptBatches, PARAM = NoCorrectParam())
        altExp(geneCorrected, "transcript") <- transcriptMerged
    }

    geneCorrected <- record_experiment_data(geneCorrected, experiment_name = "integrated", organism = organism)

    return(geneCorrected)
}

#' Reintegrate (filtered) SingleCellExperiment objects
#'
#' This function takes a SCE object and performs the below steps
#' 1) split by batch
#' 2) integrate
#' 3) run integration pipeline and save
#'
#' @param object A SingleCellExperiment objects
#' @param suffix to be appended to file saved in output dir
#' @param reduction to use default is pca
#' @param ... extra args passed to object_integration_pipeline
#'
#' @return a SingleCellExperiment object
reintegrate_object <- function(object, suffix = "", reduction = "PCA", ...) {
    organism <- metadata(object)$experiment$organism
    experiment_name <- metadata(object)$experiment$experiment_name
    objects <- splitByCol(object, "batch")
    object <- object_integration_pipeline(objects, suffix = suffix, ...)
    object <- record_experiment_data(object, experiment_name, organism)
}
