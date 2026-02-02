#' Creating Input Data for MetaHD When Individual-Level Data are Available
#'
#' @description
#' The MetaHDInput function creates input data \code{Y} (treatment effects) and \code{Slist} (within-study covariance matrices) for \code{MetaHD} when individual-level data are available. Assuming that the individual-level data are in the following format, with 'study' in column 1, 'group' in column 2 and outcomes in rest of the columns, with samples in rows.
#'
#' @usage MetaHDInput(data)
#'
#' @param data a dataframe consisting of individual-level data in the format, where 'study' in column 1, 'group' in column 2 and outcomes in rest of the columns and samples in rows.
#'
#' @return A list of objects containing :
#' \itemize{
#'   \item \code{Y}: A \eqn{K \times N} matrix of treatment effect sizes, where \eqn{K} is the number of studies and \eqn{N} is the number of outcomes.
#'   \item \code{Slist}: A list of length \eqn{K} containing \eqn{N \times N} within-study variance–covariance matrices of the treatment effects.
#' }
#'
#' @examples
#' # CREATE INPUT DATA
#' input_data <- MetaHDInput(realdata)
#'
#' ## treatment effect-sizes
#' Y <- input_data$Y
#' head(Y)
#'
#' ## within-study variance–covariance matrices
#' Slist <- input_data$Slist
#' head(Slist[[1]])
#'
#' @export MetaHDInput
#'
#' @importFrom dplyr %>% group_by summarise across everything arrange desc all_of
#' @importFrom metafor escalc
#' @importFrom tidyr gather
#' @importFrom Matrix nearPD
#' @importFrom matrixcalc is.positive.definite


MetaHDInput <- function(data){
  data <- as.data.frame(data)
  if (!is.factor(data[, 1]) || !is.factor(data[, 2])) {
    stop("Require study and group names as factors in the first and second columns respectively.")
  }
  if(length(unique(data[,1])) < 2){
    stop("Require at least two studies to prepare input data for the meta-analysis.\nEnsure that the first column contains the study names, the second column contains the groups.")
  }
  if(length(unique(data[,2]))!=2){
    stop("Restrict to two groups only.\nEnsure that the first column contains the study names, the second column contains the groups.")
  }
  if (any(is.na(data[,-c(1,2)]))) {
    stop("The dataset contains missing values. MetaHDInput requires a complete data matrix.")
  }
  names(data)[1:2] <- c("study", "group")
  group <- unique(data$group)
  N <- ncol(data[-c(1,2)])
  K <- length(unique(data[,1]))
  var_names <- names(data[-c(1,2)])
  split_data <- split(data,data$study)
  sum_data <- data %>% group_by(study, group) %>%
    summarise(across(everything(), list(Mean = ~mean(.), Sd = ~sd(.), N = ~length(.)), .names = "{fn}_{col}"),.groups = "drop") %>%
    arrange(desc(group))
  stat_data <- as.data.frame(sum_data[-c(1,2)])
  study <- unique(sum_data$study)
  meta.data <- list()
  effects <- list()
  variances <- list()
  for (i in 1:N) {
    mean_col <- (i - 1) * 3 + 1
    sd_col <- mean_col + 1
    n_col <- mean_col + 2
    meta.data[[i]] <- escalc(measure = "ROM",
                             m1i = stat_data[1:K, mean_col],
                             m2i = stat_data[(K+1):(2*K), mean_col],
                             n1i = stat_data[1:K, n_col],
                             n2i = stat_data[(K+1):(2*K), n_col],
                             sd1i = stat_data[1:K, sd_col],
                             sd2i = stat_data[(K+1):(2*K), sd_col],
                             append = FALSE)
    effects[[i]] <- meta.data[[i]][1]
    variances[[i]] <- meta.data[[i]][2]
  }
  Effects <- as.data.frame(effects)
  Variances <- as.data.frame(variances)
  colnames(Effects) <- var_names
  rownames(Effects) <- study
  colnames(Variances) <- var_names
  rownames(Variances) <- study
  var_df <- Variances
  var_df$study <- study
  var_df_long <- gather(var_df, key = "outcome", value = "var_est", all_of(var_names), factor_key=TRUE)
  sd_split <- split(sqrt(var_df_long$var_est),var_df_long$study)
  Sk <- list()
  wscormat.shrink <- list()
  for (k in 1:K) {
    wscormat.shrink[[k]] <- estimateCorMat(log(split_data[[k]][,3:(N+2)]))
    Sk[[k]] <- getCovMat(sd_split[[k]],wscormat.shrink[[k]])
    rownames(Sk[[k]]) <- colnames(Sk[[k]]) <- var_names
    if (!is.positive.definite(Sk[[k]])) {
      Sk[[k]] <- as.matrix(nearPD(Sk[[k]],keepDiag = TRUE)$mat)
    }
  }
  return(list(Y = as.matrix(Effects),
              Slist = Sk))
}

