#' Calculate Statistics
#' 
#' @description Calculates summary statistics for a numeric vector.
#' 
#' @param x Numeric vector
#' @param na.rm Whether to remove NA values (default: TRUE)
#' 
#' @return List of summary statistics
#' 
#' @export
calculate_stats <- function(x, na.rm = TRUE) {
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  
  list(
    mean = mean(x, na.rm = na.rm),
    median = median(x, na.rm = na.rm),
    sd = sd(x, na.rm = na.rm),
    min = min(x, na.rm = na.rm),
    max = max(x, na.rm = na.rm),
    n = length(x),
    n_na = sum(is.na(x))
  )
}

#' Perform T-Test
#' 
#' @description Performs a t-test between two groups.
#' 
#' @param x Numeric vector for first group
#' @param y Numeric vector for second group
#' @param paired Whether to perform a paired t-test (default: FALSE)
#' @param alternative Alternative hypothesis (default: "two.sided")
#' 
#' @return List of t-test results
#' 
#' @export
perform_ttest <- function(x, y, paired = FALSE, alternative = "two.sided") {
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric")
  }
  
  test_result <- t.test(x, y, paired = paired, alternative = alternative)
  
  list(
    statistic = test_result$statistic,
    p_value = test_result$p.value,
    alternative = test_result$alternative,
    method = test_result$method,
    estimate = test_result$estimate
  )
}

#' Calculate Correlation
#' 
#' @description Calculates correlation between two numeric vectors.
#' 
#' @param x Numeric vector
#' @param y Numeric vector
#' @param method Correlation method (default: "pearson")
#' 
#' @return List of correlation results
#' 
#' @export
calculate_correlation <- function(x, y, method = "pearson") {
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric")
  }
  
  cor_result <- cor.test(x, y, method = method)
  
  list(
    estimate = cor_result$estimate,
    p_value = cor_result$p.value,
    method = cor_result$method,
    alternative = cor_result$alternative
  )
}

#' Perform ANOVA
#' 
#' @description Performs one-way ANOVA.
#' 
#' @param formula Formula specifying the model
#' @param data Data frame containing the variables
#' 
#' @return List of ANOVA results
#' 
#' @export
perform_anova <- function(formula, data) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  
  anova_result <- aov(formula, data = data)
  summary_result <- summary(anova_result)
  
  list(
    f_statistic = summary_result[[1]]$`F value`,
    p_value = summary_result[[1]]$`Pr(>F)`,
    df = summary_result[[1]]$Df,
    residuals = summary_result[[1]]$`Mean Sq`
  )
} 