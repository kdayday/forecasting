

#' Load data from a NETCDF file ensemble forecasts
#' Assumed file dimensions: Day x Hour x Site x Lead time x Ensemble member
#' Returns an array of data: [site x member x time points]
#' Site selection and member selection selection can all be vectors of non-consecutive values
#' Time-point selection is a consecutive sequence
#' @param fname file name
#' @param members A vector of member indices
#' @param sites A vector of sites
#' @param lead_time Forecast lead time
#' @param date_start A lubridate: Start date of data to load
#' @param date_end A lubridate: End date of data to load
#' @param truncate Boolean: Whether or not to truncate the forecasts at the site maximum power
#' @param site_max_power A vector of the maximum power at ALL sites (not just those listed in sites)
get_netcdf_forecast_data <- function(fname, members, sites, lead_time, date_start, date_end,
                                     date_data_start=lubridate::ymd(20160101),
                                     ts_per_day=24, vname="power",
                                     truncate=F, site_max_power=NA) {

  if (truncate & all(is.na(site_max_power))) stop("Vector of site maximum powers required to truncate forecasts.")

  # Calculate netcdf date constants
  ndays <- interval(date_start, date_end)/days(1) + 1
  start_day <- interval(date_data_start, date_start)/days(1) + 1

  dim_counts <- c(ndays, ts_per_day, 1, 1, 1)

  # Open file
  nc <- nc_open(fname)

  # Get a matrix for this member and site
  get_member_data <- function(member, site) {
    dim_starts <- c(start_day,1, site, lead_time, member)
    return(ncvar_get(nc, varid=vname, start=dim_starts, count=dim_counts))
  }

  # Get a [member x day x hour] matrix at this site
  get_site_data <- function(site) {
    m <- sapply(members, FUN = get_member_data, site=site, simplify ="array")
    m[which(m > site_max_power[site])] <- site_max_power[site]
    return(m)
    }

  data <- sapply(sites, get_site_data, simplify="array")

  # Close the file!
  nc_close(nc)

  # Return in dimensions [site x member x day x hour]
  return(aperm(data, c(4, 3, 1, 2)))
}



