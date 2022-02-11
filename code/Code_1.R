###############
### Setup ####
##############

library(fmxdat)
library(devtools)
library(tidyverse)
library(dplyr)
library(tbl2xts)
library(rmsfuns)
library(RiskPortfolios)
library(factoextra)
library(FactoMineR)
library(lubridate)
library(MTS)
library(vars)
library(gt)
library(glue)

T40 <- readRDS("./data/T40.rds")

###################
### Wrangle 1 ####
##################

Top20_tickers <- c("ABG", "AGL", "ANG", "APN", "BHP", "CFR", "FSR", "GRT", "INL", "INP", "MNP", "MTN", "NED", "NPN", "REM",
           "SBK", "SHP", "SLM", "SOL", "VOD") # Selected based on manual inspection.

######### Equally-weighted Top 20 -> 2013.01.01 - 2021.10.29

Top_20 <- T40 %>%
dplyr::select(-J400, -J200, -Index_Name, -Short.Name, -Sector,) %>%
mutate(Tickers =gsub(" SJ Equity", "", Tickers)) %>%
filter(date > ymd(20121212)) %>%
filter(Tickers %in% Top20_tickers)

#################
### Figure 1 ####
################

facet_plot <- Top_20 %>% ggplot() + geom_line(aes(date, Return)) + facet_wrap(~Tickers) + theme_bw() +
    fmxdat::fmx_cols() +
    fmxdat::theme_fmx(title.size = ggpts(50), CustomCaption = F) +
    labs(x = "Date", y = "Returns", title = "Return series of each stock in portfolio", cex.axis = 0.5)
facet_plot


#############################
### Figure 2  & Wrangle ####
############################

Top_20_centered <- Top_20 %>%
mutate(MReturn = Return - mean(Return)) %>%
dplyr::select(-Return) %>%
rename(Return = MReturn)

tail(Top_20) # last date is 2021.10.29

return_mat_Top20 <- Top_20_centered %>% spread(Tickers, Return)

indx <- apply(return_mat_Top20, 2, function(x) any(is.na(x) | is.infinite(x))) #check for Nas.
indx # there are some NAs in ANG & APN, let's impute these.

impute_missing_returns <- function(return_mat, impute_returns_method = "NONE",
                                   Seed = 1234){
    # Make sure we have a date column called date:
    if( !"date" %in% colnames(return_mat) ) stop("No 'date' column
provided in return_mat. Try again please.")

    # Note my use of 'any' below...
    # Also note that I 'return' return_mat - which stops the function and returns return_mat.
    if( impute_returns_method %in% c("NONE", "None", "none") ) {
        if( any(is.na(return_mat)) ) warning("There are missing values in the return matrix
Consider maybe using impute_returns_method =
'Drawn_Distribution_Own' / 'Drawn_Distribution_Collective'")
        return(return_mat)
    }

    if( impute_returns_method  == "Average") {
        return_mat <-
            return_mat %>% gather(Stocks, Returns, -date) %>%
            group_by(date) %>%
            mutate(Avg = mean(Returns, na.rm=T)) %>%
            mutate(Avg = coalesce(Avg, 0)) %>%
            ungroup() %>%
            mutate(Returns = coalesce(Returns, Avg)) %>% select(-Avg) %>%
            spread(Stocks, Returns)

    } else

        if( impute_returns_method  == "Drawn_Distribution_Own") {

            set.seed(Seed)
            N <- nrow(return_mat)
            return_mat <- left_join(return_mat %>% gather(Stocks, Returns, -date),
                                    return_mat %>% gather(Stocks, Returns, -date) %>% group_by(Stocks) %>%
                                        do(Dens = density(.$Returns, na.rm=T)) %>%
                                        ungroup() %>% group_by(Stocks) %>% # done to avoid warning.
                                        do(Random_Draws = sample(.$Dens[[1]]$x, N, replace = TRUE, prob=.$Dens[[1]]$y)),
                                    by = "Stocks") %>%
                group_by(Stocks) %>% mutate(Row = row_number()) %>%
                mutate(Returns = coalesce(Returns, Random_Draws[[1]][Row])) %>%
                dplyr::select(-Random_Draws, -Row) %>% ungroup() %>% spread(Stocks, Returns)

        } else

            if( impute_returns_method  == "Drawn_Distribution_Collective") {
                set.seed(Seed)
                NAll <- nrow(return_mat %>% gather(Stocks, Returns, -date))

                return_mat <-
                    bind_cols(
                        return_mat %>% gather(Stocks, Returns, -date),
                        return_mat %>% gather(Stocks, Returns, -date) %>%
                            do(Dens = density(.$Returns, na.rm=T)) %>%
                            do(Random_Draws = sample(.$Dens[[1]]$x, NAll, replace = TRUE,
                                                     prob=.$Dens[[1]]$y)) %>% unnest(Random_Draws)) %>%
                    mutate(Returns = coalesce(Returns, Random_Draws)) %>%
                    dplyr::select(-Random_Draws) %>% spread(Stocks, Returns)

            } else

                if( impute_returns_method  == "Zero") {
                    warning("This is probably not the best idea but who am I to judge....")
                    return_mat[is.na(return_mat)] <- 0

                } else
                    stop("Please provide a valid impute_returns_method method.
Options include:\n'Average', 'Drawn_Distribution_Own',
'Drawn_Distribution_Collective' and 'Zero'.")
}
options(scipen = 999)
return_mat_Top_20_new <- impute_missing_returns(return_mat,
impute_returns_method = "Drawn_Distribution_Collective",
Seed = as.numeric(format( Sys.time(), "%Y%d%H%M")))

indx_2 <- apply(return_mat_Top_20_new, 2, function(x) any(is.na(x) | is.infinite(x))) #check for Nas.
indx_2 # No more NAs.

# Drop date column for this
return_mat_Nodate <- data.matrix(return_mat_Top_20_new[, -1])

pca <- prcomp(return_mat_Nodate, center = FALSE, scale. = FALSE) # not needed to center as all in same scale.
summary(pca)

get_pca_var(pca)$coord

scree_1 <- fviz_screeplot(pca, ncp = 12)
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 38), main = "Scree plot of PCAs", xlab = "PCA")
scree_1

contrib_1 <- fviz_contrib(pca, choice = "var", sort.val = "desc", top = 12, title = "Contributions of each share to PCA_1")

eig.val <- get_eigenvalue(pca)
eig.val

pca$rotation

#################
### Table 1 ####
################

eigenv <- get_pca_var(pca)$coord %>% as.data.frame() %>% as_tibble() %>%
dplyr::select("Dim.1", "Dim.2", "Dim.3", "Dim.4", "Dim.5", "Dim.6", "Dim.7", "Dim.8", "Dim.9", "Dim.10", "Dim.11", "Dim.12")
eigenv_table <- eigenv
eigenv_table$Ticker = Tickers
eigenv_table <- eigenv_table[, c(13,1:12)]

eigenv_table_final <- eigen_table %>%
gt() %>%
tab_header(title = "Eigenvectors of PCA Process") %>%
cols_width(everything() ~ px(40)) %>% tab_options( table.font.size = px(10))
eigenv_table_final

####################
### PVCA Model #####
###################

VARorder(return_mat_Nodate)
# AIC suggests 2 lags
# first test to see that no serial autocorrelation is present.
VAR1 <- VAR(return_mat_Nodate, p=2)
vars::serial.test(VAR1) # reject null of no serial correlation; this is not fixed by more lags.
residuals(VAR1)

m1=comVol(return_mat_Nodate, p = 2)
names(m1)


#################
### Figure 4 ####
################

eval <- round(m1$values,1)
sum(eval)
prop = round((eval/sum(eval)*100),1)

prop_calc <- tibble(Loadings = prop) %>% mutate(PC = paste0("PC_", row_number())) %>%
filter(Loadings >= 1)

prop_calc[, "PC"][[1]] <- factor(prop_calc[, "PC"][[1]], levels = prop_calc$PC)

prop_calc_plot <- prop_calc %>%
ggplot() + geom_bar(aes(PC, Loadings), stat = "identity", fill = "steelblue") +
geom_text(aes(PC, Loadings, label = Loadings), stat = "identity", vjust = 1.5, colour = "white") +
labs(x = "Principal Components", y = "Loadings", title = "Eigenvalue proportions")

####################
### Table 3 #######
###################

eigenvectors = (round(m1$M,3)) %>% as.data.frame() %>% as_tibble() %>%
dplyr::select("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")

Tickers <- T40 %>%
dplyr::select(-J400, -J200, -Index_Name, -Short.Name, -Sector,) %>%
mutate(Tickers =gsub(" SJ Equity", "", Tickers)) %>%
filter(date > ymd(20121212)) %>%
filter(Tickers %in% Top20_tickers) %>%
pull(Tickers) %>% unique()

eigen_table <- eigenvectors
eigen_table$Ticker = Tickers
eigen_table <- eigen_table[, c(13,1:12)]

eigen_table_final <- eigen_table %>%
gt() %>%
tab_header(title = "Eigenvectors of PVCA Process") %>%
cols_width(everything() ~ px(40)) %>% tab_options( table.font.size = px(10))
eigen_table_final

####################
### Figure 5 #######
###################

indiv_contrib <- eigen_table %>% dplyr::select(Ticker, V1) %>%
mutate(Contribution = (V1/sum(abs(V1))*100)) %>% dplyr:select(-V1)

indiv_contrib$Contribution = round(indiv_contrib$Contribution, 2)

indiv_contrib$Contribution = abs(indiv_contrib$Contribution)

indiv_contrib_plot <- indiv_contrib %>%
ggplot() + geom_bar(aes(Ticker, Contribution), stat = "identity", fill = "steelblue") +
geom_text(aes(Ticker, Contribution, label = Contribution), stat = "identity", vjust = 0, colour = "black") +
labs(x = "Tickers", y = "Contribution (%)", title = "Contribution of each stock to PVCA 1")

####################
### Figure 6 #######
###################

rt = m1$residuals%*%m1$M %>% as.data.frame() %>% as_tibble() %>% mutate(Period = 1:2215 )%>% dplyr::select(Period, V1, V5, V15, V19) %>%
gather(PVCA, Value, -Period)

rt_plot <- rt %>%
ggplot() + geom_line(aes(Period, Value)) + facet_wrap(~PVCA) + theme_bw() +
fmxdat::fmx_cols() +
fmxdat::theme_fmx(title.size = ggpts(50), CustomCaption = F) +
labs(x = "Period", y = "Returns", title = "Time series of selected PVCAs", cex.axis = 0.5)

#####################################################################




