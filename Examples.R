## Necessary librarires
#library(readxl)
#library(tidyr)
#library(Matrix)
#library(matrixcalc)
#library(corpcor)
################################################# EXAMPLE 1 ################################################################
### APPLICATION TO REAL DATA ###

# Read the data
GCLC_data <- read_excel("GC-LC_data.xlsx",sheet = "GC-LC_data")
Effects <- read_excel("GC-LC_data.xlsx", sheet = "Effect_sizes")
Variances <- read_excel("GC-LC_data.xlsx", sheet = "Associated_variances")

# Creating within-study covariance matrices
study <- unique(GCLC_data$Study)
K <- length(study)
var_df <- data.frame(Variances, study=study)
var_df_long <- gather(var_df, metabolite, var_est,
                      Cysteine:Valine,
                      factor_key=TRUE)
sd_split <- split(sqrt(var_df_long$var_est),var_df_long$study)

Sk <- list()
wscormat.shrink <- list()
split_data <- split(GCLC_data,GCLC_data$Study)
for (k in 1:K) {
  # within-study correlation from real data
  wscormat.shrink[[k]] <- estimateCorMat(log(split_data[[k]][,3:16]))
  # within-study covariance matrices
  Sk[[k]] <- getCovMat(sd_split[[k]],wscormat.shrink[[k]])
  if (!is.positive.definite(Sk[[k]])) {
    Sk[[k]] <- as.matrix(nearPD(Sk[[k]],keepDiag = TRUE)$mat)
  }
}

Y <- as.matrix(Effects)
Slist <- Sk

## MULTIVARIATE RANDOM-EFFECTS META-ANALYSIS, ESTIMATED WITH REML
model <- MetaHD(Y, Slist)
model$estimate

################################################# EXAMPLE 2 ################################################################
### APPLICATION TO SIMULATED DATA ###

## COMPLETE DATA EXAMPLE ##

# Read the data
Effects <- read_excel("simulated_data_1.xlsx", sheet = "effect_sizes")
WSVar <- read_excel("simulated_data_1.xlsx", sheet = "within-study_variances")
WSCorr <- read_excel("simulated_data_1.xlsx", sheet = "within-study_correlations") # UPPER TRIANGULAR ELEMENTS
K <- nrow(Effects)
var_names <- colnames(Effects)

# Creating within-study covariance matrices
Sk_var_df_initial <- data.frame(WSVar,study=1:K)
Sk_var_df_long_initial <- gather(Sk_var_df_initial, metabolite, var_est,
                                 all_of(var_names),
                                 factor_key=TRUE)
Sk_sd_initial<- split(sqrt(Sk_var_df_long_initial$var_est),Sk_var_df_long_initial$study)
Sk <- list()
for (k in 1:K) {
  Sk[[k]] <- getCovMat(Sk_sd_initial[[k]],WSCorr[k,])
  if (!is.positive.definite(Sk[[k]])) {
    Sk[[k]] <- as.matrix(nearPD(Sk[[k]],keepDiag = TRUE,maxit = 500)$mat)
  }
}

Y <- as.matrix(Effects)
Slist <- Sk

## MULTIVARIATE RANDOM-EFFECTS META-ANALYSIS FOR SIMULATED DATA (HIGH-WITHIN STUDY CORRELATIONS), ESTIMATED WITH REML
model_1 <- MetaHD(Y,
                  Slist,
                  method = "reml",
                  bscov = "unstructured")
model_1$estimate # TRUE VALUES ARE: 5, 1, 3, 2, 4, 2, 1, 1, 5, 1, 1, 5, 5, 2, 3, 5, 4, 5, 1, 1, 2, 5, 1, 2, 1, 5, 4, 1, 4, 2

################################################# EXAMPLE 3 ################################################################
## IN THE PRESENCE OF MISSING VALUES ##

# Read the data
Effects <- read_excel("simulated_data_2.xlsx", sheet = "effect_sizes")
WSVar <- read_excel("simulated_data_2.xlsx", sheet = "within-study_variances")
WSCorr <- read_excel("simulated_data_2.xlsx", sheet = "within-study_correlations") # UPPER TRIANGULAR ELEMENTS
K <- nrow(Effects)
N <- ncol(Effects)
var_names <- colnames(Effects)

# Creating within-study covariance matrices
Sk_var_df_initial <- data.frame(WSVar,study=1:K)
Sk_var_df_long_initial <- gather(Sk_var_df_initial, metabolite, var_est,
                                 all_of(var_names),
                                 factor_key=TRUE)
Sk_sd_initial<- split(sqrt(Sk_var_df_long_initial$var_est),Sk_var_df_long_initial$study)
Sk <- list()
for (k in 1:K) {
  Sk[[k]] <- getCovMat(Sk_sd_initial[[k]],WSCorr[k,])
}

Y <- as.matrix(Effects)
Slist <- Sk

## MULTIVARIATE RANDOM-EFFECTS META-ANALYSIS (HIGH-WITHIN STUDY CORRELATIONS), ESTIMATED WITH REML
model_2 <- MetaHD(Y,
                  Slist,
                  method = "reml",
                  bscov = "unstructured",
                  impute.na = TRUE)
model_2$estimate # TRUE VALUES ARE: 5, 1, 3, 2, 4, 2, 1, 1, 5, 1, 1, 5, 5, 2, 3, 5, 4, 5, 1, 1, 2, 5, 1, 2, 1, 5, 4, 1, 4, 2

## UNIVARIATE RANDOM-EFFECTS META-ANALYSIS
est <- c()
for (i in 1:N) {
  model_3 <- MetaHD(Y = na.omit(Y[,i]),
                    Slist = unlist(apply(na.omit(as.matrix(WSVar[,i])), MARGIN = 1, FUN = list), recursive = FALSE),
                    method = "reml",
                    bscov = "diag")
  est[i] <- model_3$estimate
}
est # TRUE VALUES ARE: 5, 1, 3, 2, 4, 2, 1, 1, 5, 1, 1, 5, 5, 2, 3, 5, 4, 5, 1, 1, 2, 5, 1, 2, 1, 5, 4, 1, 4, 2
