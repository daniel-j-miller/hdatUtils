

#' add_prop_CI)cold
#'
#' @param .data dataframe to add columns to
#' @param numerator the name of the numerator column
#' @param denominator the name of the denominator column
#' @param conf.level confidence level for CIs
#' @param method type of CI to be calculated
#' @param rand random seed
#' @param R output format
#'
#' @returns three dataframe columns

#' @importFrom dplyr mutate
#' @importFrom purrr map2
#' @export
prop_ci_cols <- function (.data, numerator = "n", denominator = "d", conf.level = 0.95, method = "wald", rand = 123,
                          R = 9999, bootci.type = "all", alternative = c("two.sided",
                                                                         "less", "greater"), ...)
{
  if (!is.na(pmatch(method, "wilson")))
    method <- "wilson"
  METHODS <- c("wald", "wilson", "agresti-coull", "jeffreys",
               "modified wilson", "modified jeffreys", "clopper-pearson",
               "arcsine", "logit", "witting", "wald-cc", "boot")
  method <- pmatch(method, METHODS)
  alternative <- match.arg(alternative)
  if (is.na(method))
    stop("invalid method")
  if (method == -1)
    stop("ambiguous method")
  if (length(conf.level) != 1)
    stop("'conf.level' has to be of length 1 (confidence level)")
  if (conf.level < 0.5 | conf.level > 1)
    stop("'conf.level' has to be in [0.5, 1]")
  stopifnot(R >= 1)
  R <- trunc(R)
  alpha <- 1 - conf.level

  dplyr::mutate(.data, purrr::map2(.data[[numerator]], .data[[denominator]], function(n,d) {

    numerator <- as.integer(n)
    denominator <- as.integer(d)

    if (alternative != "two.sided")
      alpha <- 2 * alpha
    kappa <- qnorm(1 - alpha/2)
    p.hat <- numerator/denominator
    q.hat <- 1 - p.hat
    Infos <- NULL
    if (method == 1) {
      est <- p.hat
      term2 <- kappa * sqrt(p.hat * q.hat)/sqrt(denominator)
      CI.lower <- max(0, p.hat - term2)
      CI.upper <- min(1, p.hat + term2)
      Infos <- term2/kappa
      names(Infos) <- "standard error of prob"
    }
    if (method == 2) {
      est <- p.hat
      term1 <- (numerator + kappa^2/2)/(denominator + kappa^2)
      term2 <- kappa * sqrt(denominator)/(denominator + kappa^2) * sqrt(p.hat *
                                                                          q.hat + kappa^2/(4 * denominator))
      CI.lower <- max(0, term1 - term2)
      CI.upper <- min(1, term1 + term2)
      Infos <- term2/kappa
      names(Infos) <- "standard error of prob"
    }
    if (method == 3) {
      x.tilde <- numerator + kappa^2/2
      n.tilde <- denominator + kappa^2
      p.tilde <- x.tilde/n.tilde
      q.tilde <- 1 - p.tilde
      est <- p.tilde
      term2 <- kappa * sqrt(p.tilde * q.tilde)/sqrt(n.tilde)
      CI.lower <- max(0, p.tilde - term2)
      CI.upper <- min(1, p.tilde + term2)
      Infos <- term2/kappa
      names(Infos) <- "standard error of prob"
    }
    if (method == 4) {
      est <- p.hat
      if (numerator == 0)
        CI.lower <- 0
      else CI.lower <- qbeta(alpha/2, numerator + 0.5, denominator - numerator + 0.5)
      if (numerator == denominator)
        CI.upper <- 1
      else CI.upper <- qbeta(1 - alpha/2, numerator + 0.5, denominator - numerator +
                               0.5)
    }
    if (method == 5) {
      est <- p.hat
      term1 <- (numerator + kappa^2/2)/(denominator + kappa^2)
      term2 <- kappa * sqrt(denominator)/(denominator + kappa^2) * sqrt(p.hat *
                                                                          q.hat + kappa^2/(4 * denominator))
      if ((denominator <= 50 & numerator %in% c(1, 2)) | (denominator >= 51 & numerator %in% c(1:3)))
        CI.lower <- 0.5 * qchisq(alpha, 2 * numerator)/denominator
      else CI.lower <- max(0, term1 - term2)
      if ((denominator <= 50 & numerator %in% c(denominator - 1, denominator - 2)) | (denominator >= 51 & numerator %in%
                                                                                      c(denominator - (1:3))))
        CI.upper <- 1 - 0.5 * qchisq(alpha, 2 * (denominator - numerator))/denominator
      else CI.upper <- min(1, term1 + term2)
      Infos <- term2/kappa
      names(Infos) <- "standard error of prob"
    }
    if (method == 6) {
      est <- p.hat
      if (numerator == denominator)
        CI.lower <- (alpha/2)^(1/denominator)
      else {
        if (numerator <= 1)
          CI.lower <- 0
        else CI.lower <- qbeta(alpha/2, numerator + 0.5, denominator - numerator +
                                 0.5)
      }
      if (numerator == 0)
        CI.upper <- 1 - (alpha/2)^(1/denominator)
      else {
        if (numerator >= denominator - 1)
          CI.upper <- 1
        else CI.upper <- qbeta(1 - alpha/2, numerator + 0.5, denominator -
                                 numerator + 0.5)
      }
    }
    if (method == 7) {
      est <- p.hat
      CI.lower <- qbeta(alpha/2, numerator, denominator - numerator + 1)
      CI.upper <- qbeta(1 - alpha/2, numerator + 1, denominator - numerator)
    }
    if (method == 8) {
      p.tilde <- (numerator + 0.375)/(denominator + 0.75)
      est <- p.tilde
      CI.lower <- sin(asin(sqrt(p.tilde)) - 0.5 * kappa/sqrt(denominator))^2
      CI.upper <- sin(asin(sqrt(p.tilde)) + 0.5 * kappa/sqrt(denominator))^2
    }
    if (method == 9) {
      est <- p.hat
      lambda.hat <- log(numerator/(denominator - numerator))
      V.hat <- denominator/(numerator * (denominator - numerator))
      lambda.lower <- lambda.hat - kappa * sqrt(V.hat)
      lambda.upper <- lambda.hat + kappa * sqrt(V.hat)
      CI.lower <- exp(lambda.lower)/(1 + exp(lambda.lower))
      CI.upper <- exp(lambda.upper)/(1 + exp(lambda.upper))
    }
    if (method == 10) {
      set.seed(rand)
      x.tilde <- numerator + runif(1, min = 0, max = 1)
      pbinom.abscont <- function(q, size, prob) {
        v <- trunc(q)
        term1 <- pbinom(v - 1, size = size, prob = prob)
        term2 <- (q - v) * dbinom(v, size = size, prob = prob)
        return(term1 + term2)
      }
      qbinom.abscont <- function(p, size, x) {
        fun <- function(prob, size, x, p) {
          pbinom.abscont(x, size, prob) - p
        }
        uniroot(fun, interval = c(0, 1), size = size, x = x,
                p = p)$root
      }
      est <- p.hat
      CI.lower <- qbinom.abscont(1 - alpha, size = denominator, x = x.tilde)
      CI.upper <- qbinom.abscont(alpha, size = denominator, x = x.tilde)
    }
    if (method == 11) {
      est <- p.hat
      term2 <- kappa * sqrt(p.hat * q.hat)/sqrt(denominator)
      CC <- 0.5/denominator
      CI.lower <- max(0, p.hat - term2 - CC)
      CI.upper <- min(1, p.hat + term2 + CC)
      Infos <- term2/kappa
      names(Infos) <- "standard error of prob"
    }
    if (method == 12) {
      if (numerator == 0 | numerator == denominator)
        warning("All observations are identical.\n", "Choose a different method for computing the confidence interval!")
      est <- p.hat
      DATA <- numeric(denominator)
      DATA[1:numerator] <- 1
      boot.rf <- function(x, i) {
        p <- mean(x[i])
        n <- length(i)
        c(p, p * (1 - p)/n)
      }
      boot.out <- boot(DATA, statistic = boot.rf, R = R, ...)
      CI <- try(boot.ci(boot.out, type = bootci.type, conf = 1 -
                          alpha), silent = TRUE)
      if (inherits(CI, "try-error"))
        stop("Function 'boot.ci' returned an error. Please try a different 'bootci.type'.")
      Infos <- c(sqrt(p.hat * q.hat)/sqrt(denominator), sqrt(var(boot.out$t[,
                                                                            1])))
      names(Infos) <- c("standard error of prob", "bootstrap standard error of prob")
    }
    if (alternative == "less")
      if (method == 12) {
        if ("normal" %in% names(CI)) {
          CI$normal[1, 1] <- conf.level
          CI$normal[1, 2] <- 0
        }
        if ("basic" %in% names(CI)) {
          CI$basic[1, 1] <- conf.level
          CI$basic[1, 4] <- 0
        }
        if ("student" %in% names(CI)) {
          CI$student[1, 1] <- conf.level
          CI$student[1, 4] <- 0
        }
        if ("percent" %in% names(CI)) {
          CI$percent[1, 1] <- conf.level
          CI$percent[1, 4] <- 0
        }
        if ("bca" %in% names(CI)) {
          CI$bca[1, 1] <- conf.level
          CI$bca[1, 4] <- 0
        }
      }
    else CI.lower <- 0
    if (alternative == "greater")
      if (method == 12) {
        if ("normal" %in% names(CI)) {
          CI$normal[1, 1] <- conf.level
          CI$normal[1, 3] <- 1
        }
        if ("basic" %in% names(CI)) {
          CI$basic[1, 1] <- conf.level
          CI$basic[1, 5] <- 1
        }
        if ("student" %in% names(CI)) {
          CI$student[1, 1] <- conf.level
          CI$student[1, 5] <- 1
        }
        if ("percent" %in% names(CI)) {
          CI$percent[1, 1] <- conf.level
          CI$percent[1, 5] <- 1
        }
        if ("bca" %in% names(CI)) {
          CI$bca[1, 1] <- conf.level
          CI$bca[1, 5] <- 1
        }
      }
    else CI.upper <- 1
    if (method != 12) {
      stat_cols <- data.frame(p = est * 100,
                              CI_lower = CI.lower * 100,
                              CI_upper = CI.upper * 100)}}) |>
      list_rbind())

}
