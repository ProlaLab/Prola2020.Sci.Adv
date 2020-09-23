# Those function are used to transform the coefs and variance-covariance matrices of a triexponential
# model fit from reparametrized parameters (`thetas`, `Bs`) into "natural" parameters (characteristic
# times `taus`, proportions `props`).

vcovsMeans <- function (x, ...) {
  UseMethod("vcovsMeans", x)
}

vcovsMeans.default <- function(x, ...) {
  setNames(lapply(
    names(x),
    function (nm) {
      if (is.null(x[[nm]])) {
        # Did not converge
        return (NULL)
      }
      rand <- mvtnorm::rmvnorm(n=10000, mean=coef(x[[nm]]),
        sigma=vcov(x[[nm]]), method="svd")
      prop_sum <- rand[,"B1"] + rand[,"B2"] + rand[,"B3"]
      transformed_rand <- cbind(
        tau1 = exp(-rand[,"theta1"]),
        tau2 = exp(-rand[,"theta2"]),
        tau3 = exp(-rand[,"theta3"]),
        prop1 = rand[,"B1"] / prop_sum,
        prop2 = rand[,"B2"] / prop_sum,
        prop3 = rand[,"B3"] / prop_sum
      )
      list(
        cov = cov(transformed_rand),
        mean = colMeans(transformed_rand)
      )
    }
  ), names(x))
}

group.varsMeans <- function (meta, triexp.vcovsMeans) {
  triexp.varsMeans.grouped <- list()
  for (type in levels(meta$Type)) {
    for (pop in levels(meta$Population)) {
      universe <- triexp.vcovsMeans[meta$Type == type & meta$Population == pop]
      universe <- universe[!sapply(universe, is.null)]

      # Get proper ordering of parameters
      ordered_universe <- lapply(
        universe,
        function (x) {
          tau_ord <- order(x$mean[c("tau1", "tau2", "tau3")])
          all_ord <- c(tau_ord, tau_ord + 3)
          ans_cov <- x$cov[all_ord, all_ord]
          nms <- names(x$mean)
          colnames(ans_cov) <- nms
          rownames(ans_cov) <- nms
          ans_mean <- x$mean[all_ord]
          names(ans_mean) <- nms
          list(
            cov = ans_cov,
            mean = ans_mean
          )
        }
      )

      # Get mean variance
      mean_var <- rowMeans(sapply(ordered_universe, function (x) diag(x$cov)))

      # Get variance of the mean
      var_mean <- diag(cov(t(sapply(ordered_universe, function (x) x$mean))))

      # Combine
      triexp.varsMeans.grouped[[paste0(type, ".", pop)]] <- list(
        var = (mean_var + var_mean) / length(ordered_universe),
        mean = rowMeans(sapply(ordered_universe, function (x) x$mean))
      )
    }
  }

  triexp.varsMeans.grouped
}

freqPlotCov <- function (triexp.vcovsMeans) {
  df_freq_plot2 <- do.call(rbind, lapply(
    names(triexp.vcovsMeans),
    function (nm) {
      if (is.null(triexp.vcovsMeans[[nm]])) {
        return(data.frame())
      }
      alpha_level <- 0.05
      q <- qnorm(alpha_level / 2)
      ivals <- triexp.vcovsMeans[[nm]]$mean + outer(q * sqrt(diag(triexp.vcovsMeans[[nm]]$cov)), c(-1, 0, 1))
      data.frame(
        Model = with(subset(meta, File == nm), interaction(Population, Type)),
        tau_lower = ivals[c("tau1", "tau2", "tau3"), 1],
        tau_est. = ivals[c("tau1", "tau2", "tau3"), 2],
        tau_upper = ivals[c("tau1", "tau2", "tau3"), 3],
        prop_lower = ivals[c("prop1", "prop2", "prop3"), 1],
        prop_est. = ivals[c("prop1", "prop2", "prop3"), 2],
        prop_upper = ivals[c("prop1", "prop2", "prop3"), 3]
      )
    }
  ))
}

freqPlotVar <- function (triexp.varsMeans.grouped) {
  do.call(rbind, lapply(
    names(triexp.varsMeans.grouped),
    function (nm) {
      alpha_level <- 0.05
      q <- qnorm(alpha_level / 2)
      ivals <- triexp.varsMeans.grouped[[nm]]$mean + outer(q * sqrt(triexp.varsMeans.grouped[[nm]]$var), c(-1, 0, 1))
      data.frame(
        Model = nm,
        tau_lower = ivals[c("tau1", "tau2", "tau3"), 1],
        tau_est. = ivals[c("tau1", "tau2", "tau3"), 2],
        tau_upper = ivals[c("tau1", "tau2", "tau3"), 3],
        prop_lower = ivals[c("prop1", "prop2", "prop3"), 1],
        prop_est. = ivals[c("prop1", "prop2", "prop3"), 2],
        prop_upper = ivals[c("prop1", "prop2", "prop3"), 3]
      )
    }
  ))
}

plotFreq <- function (df_freq_plot2) {
  ggplot(df_freq_plot2, aes(x = tau_est., y = prop_est., color = Model))+
    geom_errorbar(aes(ymin = prop_lower, ymax = prop_upper))+
    geom_errorbarh(aes(xmin = tau_lower, xmax = tau_upper))+
    theme_classic()+
    theme(legend.position="bottom")+
    guides(color=guide_legend(nrow=3, byrow=TRUE))
}
