#' @format This package contains two data sets. The first contains NOT DONE information:
#'            var1: Description
#'            var2: Description
#'            var3: Description
#'            var4: Description
#'            var5: Description
#'
#'          The second data set contains additional NOT DONE information:
#'            var1: Description
#'            var2: Description
#'            var3: Description
#'            var4: Description
#'            var5: Description
#' @source loan (Dataset 1): NOT DONE
#' @source loan.aux (Dataset 2): NOT DONE

# Tidying loan data - NOT DONE
loan = haven::read_dta("ACCES_credit_WD_data.dta")

# Tidying loan.aux data - NOT DONE
loan.aux = haven::read_dta("puntos_corte_2000-2014.dta")

# NOT DONE: Once wrangling is finished, run this line to create new 'data' folder with clean datasets:
usethis::use_data(loan, loan.aux, internal = TRUE, overwrite = TRUE)

