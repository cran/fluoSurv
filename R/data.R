#' An injection experimental setup
#'
#' A dataset containing the description of an experimental setup where
#' 96 larvae of \emph{Galleria mellonella} have been injected with a culture
#' of the bacterium \emph{Xenorhabdus nematophila}.
#'
#' @name setup
#' @format A data frame with 96 rows and 4 variables:
#' \describe{
#'   \item{well}{Well name.}
#'   \item{dilution}{Dilution factor (log-transformed) of the injected culture. 1 therefore means 10 fold dilution, while LB corresponds to negative control where insects have been injected with sterile LB culture medium.}
#'   \item{time_injection}{Time of the injection.}
#'   \item{dead}{Is the insect dead at the end of the experiment?}
#' }
#' @format csv
NULL

#' Fluorescence data measured on \emph{Galleria mellonella}
#'
#' A dataset containing fluoresence measurements produced by a BioTek microplate
#' reader. Measurements are taken from 96 larvae of \emph{Galleria mellonella} which
#' have been injected with a culture of the bacterium \emph{Xenorhabdus nematophila}.
#'
#' @name galleria
#' @format A data frame with 382944 rows and 8 variables:
#' \describe{
#'   \item{well}{Well name.}
#'   \item{value}{Intensity measurement.}
#'   \item{t}{Time.}
#'   \item{num}{Is the insect dead at the end of the experiment?}
#'   \item{read}{Read number, usually 1. This number will be greater than one when a combination of excitation and emission wavelengths is measured several times, with different gains.}
#'   \item{exc}{Excitation wavelength.}
#'   \item{em}{Emission wavelength.}
#'   \item{ID_read}{ID of the read. Combines read number and wavelengths to produce a unique ID for each set of fluorescence measurement. For example, if GFP fluoerscence has been measured with two different gains, the two measurements will be 1_485_535 and 2_485_535.}
#' }
#' @format csv
NULL
