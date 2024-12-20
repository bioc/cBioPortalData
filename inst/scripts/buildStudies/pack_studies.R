library(cBioPortalData)

## PACK BUILD
message("PACK BUILD")

cbioportal <- cBioPortal()
studies <- stats::setNames(nm = getStudies(cbioportal)[["studyId"]])

comp_pack <- vector("logical", length(studies))
names(comp_pack) <- studies

err_pack <- vector("character", length(studies))
names(err_pack) <- studies

if (identical(tolower(Sys.getenv("IS_BIOC_BUILD_MACHINE")), "true")) {

    for (pack_stud in studies) {
        message("Working on: ", pack_stud)
        dats <- tryCatch({
            cBioDataPack(cancer_study_id = pack_stud, ask = FALSE)
        }, error = function(e) conditionMessage(e))
        success <- is(dats, "MultiAssayExperiment")
        if (success)
            comp_pack[[pack_stud]] <- success
        else
            err_pack[[pack_stud]] <- dats
        ## try to free up memory
        gc()
        ## clean up data
        removePackCache(cancer_study_id = pack_stud, dry.run = FALSE)
    }

} else if (identical(Sys.getenv("IS_SUPERMICRO_MACHINE"), "TRUE")) {

    library(BiocParallel)
    params <- MulticoreParam(
        workers = 64, stop.on.error = FALSE, progressbar = TRUE
    )

    res_pack <- bplapply(X = studies, FUN = function(x) {
        dats <- tryCatch({
            cBioPortalData::cBioDataPack(cancer_study_id = x, ask = FALSE)
        }, error = function(e) conditionMessage(e))
        comp <- is(dats, "MultiAssayExperiment")
        if (!comp)
            err <- dats
        else
            err <- ""
        list(comp_pack = comp, err_pack = err)
    }, BPPARAM = params)
    comp_pack <- vapply(res_pack, `[[`, logical(1L), "comp_pack")
    err_pack <- vapply(res_pack, `[[`, character(1L), "err_pack")

}

err_pack <- Filter(nchar, err_pack)
err_pack_info <- lapply(setNames(nm = unique(err_pack)),
    function(x) names(err_pack)[err_pack == x])
# table(err_pack)

err_pack_info_file <- "inst/extdata/pack/err_pack_info.json"
jsonlite::write_json(err_pack_info, path = err_pack_info_file)
jsonlite::fromJSON(err_pack_info_file)

pack_build <- rev(stack(comp_pack))
names(pack_build) <- c("studyId", "pack_build")
pack_build[["studyId"]] <- as.character(pack_build[["studyId"]])

pack_file_json <- "inst/extdata/pack/pack_build.json"
prev <- jsonlite::fromJSON(pack_file_json) |> as.data.frame()

if (!identical(prev, pack_build))
    pack_build |> jsonlite::write_json(path = pack_file_json, pretty = TRUE)
