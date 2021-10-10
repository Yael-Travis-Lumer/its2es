#' Monthly unemployment count and percent in Israel between January 2013 and May 2021
#'
#' A dataset containing the monthly unemployment count, the monthly labour force, and the monthly unemployment percent, as well as other attributes.
#'
#' @format A data frame with 101 rows and 8 variables:
#' \describe{
#'   \item{year}{The year, between 2013 to 2021}
#'   \item{month}{The number month, between 1 to 12}
#'   \item{unemployed}{Total monthly count of unemployed}
#'   \item{labour}{Total monthly count of labour force}
#'   \item{percent}{Monthly unemployment percent}
#'   \item{dt}{The date}
#'   \item{time}{A numeric vector counting the number of months (1 to 101)}
#' }
#' @source \url{https://www.cbs.gov.il/en/publications/Pages/2021/Labour-Force-Survey-Monthly-May-2021.aspx}
"unemployed"


#' Monthly deaths in Israel between January 2001 and May 2021
#'
#' A dataset containing the monthly number of deaths, the estimated monthly population size, and the monthly mortality percent, for all the population and for groups of different age and sex.
#'
#' @format A data frame with 245 rows and 7 variables:
#' \describe{
#'   \item{Month}{The month number, between 1 to 12}
#'   \item{Year}{The year, between 2001 to 2021}
#'   \item{monthly_total}{Total monthly count of deaths for the entire population}
#'   \item{Date}{The date}
#'   \item{monthly_est}{Monthly estimated population size}
#'   \item{percent}{Monthly mortality percent for the entire population}
#'   \item{time}{A numeric vector counting the number of months (1 to 245)}
#' }
#' @source \url{https://www.cbs.gov.il/he/Pages/search/TableMaps.aspx?CbsSubject=%D7%AA%D7%9E%D7%95%D7%AA%D7%94%20%D7%95%D7%AA%D7%95%D7%97%D7%9C%D7%AA%20%D7%97%D7%99%D7%99%D7%9D}
"Israel_mortality"


