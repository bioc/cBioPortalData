test_that("cBioPortal API is working with most studies", {

    cbio <- cBioPortal()
    studies <- getStudies(cbio)[["studyId"]]

    isMAE <- structure(vector("logical", length(studies)), .Names = studies)

    for (api_stud in studies) {
        message("Working on: ", api_stud)
        result <- try({
            cBioPortalData(
                api = cbio,
                studyId = api_stud,
                genePanelId = "IMPACT341",
                check_build = FALSE
            )
        })
        isMAE[api_stud] <- is(result, "MultiAssayExperiment")
        removeDataCache(
            api = cbio,
            studyId = api_stud,
            genePanelId = "IMPACT341",
            dry.run = FALSE
        )
    }

    successrate <- (100 * sum(isMAE)) / length(isMAE)

    expect_true(successrate > 80)
})
