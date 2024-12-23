% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{cBioCache}
\alias{cBioCache}
\alias{setCache}
\alias{removePackCache}
\title{Manage cache / download directories for study data}
\usage{
cBioCache(..., ask = interactive())

setCache(
  directory = tools::R_user_dir("cBioPortalData", "cache"),
  verbose = TRUE,
  ask = interactive()
)

removePackCache(cancer_study_id, dry.run = TRUE)
}
\arguments{
\item{...}{For \code{cBioCache}, arguments passed to \code{setCache}}

\item{ask}{logical (default TRUE when interactive session) Confirm the file
location of the cache directory}

\item{directory}{The file location where the cache is located. Once set
future downloads will go to this folder.}

\item{verbose}{Whether to print descriptive messages}

\item{cancer_study_id}{character(1) The \code{studyId} from \code{getStudies}}

\item{dry.run}{logical Whether or not to remove cache files (default TRUE).}
}
\value{
cBioCache: The path to the cache location
}
\description{
Managing data downloads is important to save disk space and
avoid re-downloading data files. This can be done via the integrated
\code{BiocFileCache} system.
}
\section{cBioCache}{

Get the directory location of the cache. It will prompt the user to create
a cache if not already created. A specific directory can be used via
\code{setCache}.
}

\section{setCache}{

Specify the directory location of the data cache. By default, it will
go to the user directory as given by:

\if{html}{\out{<div class="sourceCode">}}\preformatted{tools::R_user_dir("cBioPortalData", "cache")
}\if{html}{\out{</div>}}
}

\section{removePackCache}{

Some files may become corrupt when downloading, this function allows
the user to delete the tarball associated with a \code{cancer_study_id} in the
cache. This only works for the \code{cBioDataPack} function. To remove the entire
\code{cBioPortalData} cache, run \code{unlink("~/.cache/cBioPortalData")}.
}

\examples{

cBioCache()

removePackCache("acc_tcga", dry.run = TRUE)

}
