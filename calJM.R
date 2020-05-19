# object = jointFit3
# newdata = pbc2
# Tstart = 5
# Thoriz = 7
# idVar = "id"
# simulate = FALSE
# M = 100
# pl = TRUE
# include.cens = TRUE

calJM <- function(object, newdata, Tstart, Thoriz, idVar = "id", 
                  simulate = FALSE, M = 100, pl = TRUE, 
                  include.cens = TRUE, ...) {
  if (!inherits(object, "JMbayes")) 
    stop("Use only with 'JMbayes' objects.\n")
  if (!is.data.frame(newdata) || nrow(newdata) == 0) 
    stop("'newdata' must be a data.frame with more than one rows.\n")
  if (is.null(newdata[[idVar]])) 
    stop("'idVar' not in 'newdata.\n'")
 
  id <- newdata[[idVar]]
  id <- match(id, unique(id))
  TermsT <- object$Terms$termsT
  SurvT <- model.response(model.frame(TermsT, newdata))
  is_counting <- attr(SurvT, "type") == "counting"
  Time <- if (is_counting) {
    ave(SurvT[, 2], id, FUN = function(x) tail(x, 1))
  } else {
    SurvT[, 1]
  }
  timeVar <- object$timeVar
  newdata2 <- newdata[Time > Tstart, ]
  SurvT <- model.response(model.frame(TermsT, newdata2))
  if (is_counting) {
    id2 <- newdata2[[idVar]]
    f <- factor(id2, levels = unique(id2))
    Time <- ave(SurvT[, 2], f, FUN = function(x) tail(x, 
                                                      1))
    delta <- ave(SurvT[, 3], f, FUN = function(x) tail(x, 
                                                       1))
  } else {
    Time <- SurvT[, 1]
    delta <- SurvT[, 2]
  }
  timesInd <- newdata2[[timeVar]] <= Tstart
  aliveThoriz <- newdata2[Time > Thoriz & timesInd, ]
  deadThoriz <- newdata2[Time <= Thoriz & delta == 1 & timesInd, ]
  indCens <- Time < Thoriz & delta == 0 & timesInd
  censThoriz <- newdata2[indCens, ]
  nr <- length(unique(newdata2[[idVar]]))
  idalive <- unique(aliveThoriz[[idVar]])
  iddead <- unique(deadThoriz[[idVar]])
  idcens <- unique(censThoriz[[idVar]])

  Surv.aliveThoriz <- if (is_counting) {
    survfitJM(object, newdata = aliveThoriz, idVar = idVar, 
              simulate = simulate, M = M, survTimes = Thoriz, 
              last.time = rep(Tstart, length(idalive)), LeftTrunc_var = all.vars(TermsT)[1L])
  } else {
    survfitJM(object, newdata = aliveThoriz, idVar = idVar, 
              simulate = simulate, M = M, survTimes = Thoriz, 
              last.time = rep(Tstart, length(idalive)))
  }
  Surv.deadThoriz <- if (is_counting) {
    survfitJM(object, newdata = deadThoriz, idVar = idVar, 
              simulate = simulate, survTimes = Thoriz, last.time = rep(Tstart, 
                                                                       length(iddead)), LeftTrunc_var = all.vars(TermsT)[1L])
  } else {
    survfitJM(object, newdata = deadThoriz, idVar = idVar, 
              simulate = simulate, survTimes = Thoriz, last.time = rep(Tstart, 
                                                                       length(iddead)))
  }
  Surv.aliveThoriz <- sapply(Surv.aliveThoriz$summaries, "[", 2)
  Surv.deadThoriz <- sapply(Surv.deadThoriz$summaries, "[", 2)
  if (nrow(censThoriz)) {
    Surv.censThoriz <- if (is_counting) {
      survfitJM(object, newdata = censThoriz, idVar = idVar, 
                simulate = simulate, M = M, survTimes = Thoriz, 
                last.time = rep(Tstart, length(idcens)), LeftTrunc_var = all.vars(TermsT)[1L])
    } else {
      survfitJM(object, newdata = censThoriz, idVar = idVar, 
                simulate = simulate, M = M, survTimes = Thoriz, 
                last.time = rep(Tstart, length(idcens)))
    }
    tt <- Time[indCens]
    weights <- if (is_counting) {
      survfitJM(object, newdata = censThoriz, idVar = idVar, 
                simulate = simulate, M = M, survTimes = Thoriz, 
                last.time = tt[!duplicated(censThoriz[[idVar]])], 
                LeftTrunc_var = all.vars(TermsT)[1L])
    }
    else {
      survfitJM(object, newdata = censThoriz, idVar = idVar, 
                simulate = simulate, M = M, survTimes = Thoriz, 
                last.time = tt[!duplicated(censThoriz[[idVar]])])
    }
    Surv.censThoriz <- sapply(Surv.censThoriz$summaries, 
                              "[", 2)
    weights <- sapply(weights$summaries, "[", 2)
  } else {
    Surv.censThoriz <- weights <- NA
  }
  
  # until here I copy paste the prederrJM function from the JMbayes package

  
  statusCens <- rbinom(n = length(weights), size = 1, prob = weights)
  matCens <- data.frame(pred_prob = Surv.censThoriz, status = statusCens)
  
  # New part for calibration plot
  matDead <- data.frame(pred_prob = 1 - Surv.deadThoriz, status = 1)
  matAlive <- data.frame(pred_prob = 1 - Surv.aliveThoriz, status = 0)
  if (include.cens == TRUE) {
    mat <- rbind(matDead, matAlive, matCens)
  } else {
    mat <- rbind(matDead, matAlive)
  }


  # Logistic calibration curve
  y <- mat$status
  p <- mat$pred_prob
  f.recal <- glm(y ~ p, data = mat, family = binomial)
  logit <- seq(-7, 7, length = 200)
  prob <- plogis(logit)
  newdata = list(p = prob)
  log.cal.curve <- predict(f.recal, newdata = newdata, type = "response")
  log.cal.curve <- list(x = prob, y = log.cal.curve)
  # Harrell code - logistic calibration curve
  # y <- mat$status
  # p <- mat$pred_prob
  # logit <- qlogis(p)
  # i <- !is.infinite(logit)
  # f.recal <- lrm.fit(logit[i], y = y[i])
  # logit <- seq(-7, 7, length = 200)
  # prob <- plogis(logit)
  # pred.prob <- f.recal$coef[1] + f.recal$coef[2] * logit 
  # pred.prob2 <- plogis(pred.prob)
  
  # Harrel code - optionally a smooth nonparametric fit
  sm.nonpar.fit <- lowess(mat$pred_prob, mat$status, iter = 0)
  
  # Logistic calibration curve with splines terms
  y <- mat$status
  p <- mat$pred_prob
  f.recal <- glm(y ~ ns(p, 2), data = mat, family = binomial)
  logit <- seq(-7, 7, length = 200)
  prob <- plogis(logit)
  newdata = list(p = prob)
  log.cal.curve.poly <- predict(f.recal, newdata = newdata, type = "response")
  log.cal.curve.poly <- list(x = prob, y = log.cal.curve.poly)
  
  out <- list(log.cal.curve = log.cal.curve,
              log.cal.curve.poly = log.cal.curve.poly,
              sm.nonpar.fit = sm.nonpar.fit)
  
  if (pl == TRUE) {
    plot(0.5, 0.5, xlim = c(0, 1), ylim = c(0, 1), type = "n",
         xlab = "Predicted Probability",
         ylab = "Actual Probability",)
    abline(0, 1, lwd = 6, col = gray(0.85))
    abline(v = max(p), col = 2, lty = 2)
    lines(log.cal.curve, lty = 1, lwd = 2)
    lines(log.cal.curve.poly, lty = 2)
    lines(sm.nonpar.fit, lty = 3)
    leg <- "Ideal"
    leg <- c(leg, "Logistic calibration")
    leg <- c(leg, "Logistic calibration (splines)")
    leg <- c(leg, "Nonparametric")
    leg <- c(leg, "Max predicted prob")
    lp = c(-0.05, 1.05)
    lp = list(x = lp[1], y = lp[2])
    legend(lp, leg, lty = c(1,1,2,3,6), cex = 1, y.intersp = 0.5,
           lwd = c(6,2,1,1,1), col = c(gray(0.85), 1, 1, 1, 2), bty = "n")
  }
  out
# CHECK plot with val.prob(mat$pred_prob, mat$status)
# Still to check: polynomial logistic calibration - not sure
# Censoring is excluded for now
}


#calJM(object = jointFit1, newdata = pbc2, Tstart = 5, Thoriz = 7)
  
calJM.mv <- function(object, newdata, Tstart, Thoriz, idVar = "id", 
                  simulate = FALSE, M = 100, pl = TRUE, 
                  include.cens = TRUE, ...) {
  if (!inherits(object, "mvJMbayes"))
    stop("Use only with 'mvJMbayes' objects.\n")
  if (!is.data.frame(newdata) || nrow(newdata) == 0)
    stop("'newdata' must be a data.frame with more than one rows.\n")
  if (is.null(newdata[[idVar]]))
    stop("'idVar' not in 'newdata.\n'")
  id <- newdata[[idVar]]
  id <- match(id, unique(id))
  TermsT <- object$model_info$coxph_components$Terms
  environment(TermsT) <- parent.frame()#.GlobalEnv
  SurvT <- model.response(model.frame(TermsT, newdata)) 
  is_counting <- attr(SurvT, "type") == "counting"
  is_interval <- attr(SurvT, "type") == "interval"
  Time <- if (is_counting) {
    ave(SurvT[, 2], id, FUN = function (x) tail(x, 1))
  } else if (is_interval) {
    Time1 <- SurvT[, "time1"]
    Time2 <- SurvT[, "time2"]
    Time <- Time1
    Time[Time2 != 1] <- Time2[Time2 != 1]
    Time
  } else {
    SurvT[, 1]
  }
  timeVar <- object$model_info$timeVar
  newdata2 <- newdata[Time > Tstart, ]
  id2 <- newdata2[[idVar]]
  SurvT <- model.response(model.frame(TermsT, newdata2))
  if (is_counting) {
    f <- factor(id2, levels = unique(id2))
    Time <- ave(SurvT[, 2], f, FUN = function (x) tail(x, 1))
    event <- ave(SurvT[, 3], f, FUN = function (x) tail(x, 1))
  } else if (is_interval) {
    Time1 <- SurvT[, "time1"]
    Time2 <- SurvT[, "time2"]
    Time <- Time1
    Time[Time2 != 1] <- Time2[Time2 != 1]
    event <- SurvT[, "status"]
  } else {
    Time <- SurvT[, 1]
    event <- SurvT[, 2]
  }
  timesInd <- newdata2[[timeVar]] <= Tstart
  aliveThoriz <- newdata2[Time > Thoriz & timesInd, ]
  deadThoriz <- newdata2[Time <= Thoriz & (event == 1 | event == 3) & timesInd, ]
  indCens <- Time < Thoriz & (event == 0 | event == 2) & timesInd
  censThoriz <- newdata2[indCens, ]
  nr <- length(unique(newdata2[[idVar]]))
  idalive <- unique(aliveThoriz[[idVar]])
  iddead <- unique(deadThoriz[[idVar]])
  idcens <- unique(censThoriz[[idVar]])
  
  Surv.aliveThoriz <- if (is_counting) {
    survfitJM(object, newdata = aliveThoriz, idVar = idVar, M = M,
              survTimes = Thoriz, last.time = rep(Tstart, length(idalive)),
              LeftTrunc_var = all.vars(TermsT)[1L])
  } else {
    survfitJM(object, newdata = aliveThoriz, idVar = idVar, M = M,
              survTimes = Thoriz, last.time = rep(Tstart, length(idalive)))
    
  }
  Surv.deadThoriz <- if (is_counting) {
    survfitJM(object, newdata = deadThoriz, idVar = idVar, 
              survTimes = Thoriz, last.time = rep(Tstart, length(iddead)),
              LeftTrunc_var = all.vars(TermsT)[1L])
  } else {
    survfitJM(object, newdata = deadThoriz, idVar = idVar, 
              survTimes = Thoriz, last.time = rep(Tstart, length(iddead)))
  }
  Surv.aliveThoriz <- sapply(Surv.aliveThoriz$summaries, "[", 2)
  Surv.deadThoriz <- sapply(Surv.deadThoriz$summaries, "[", 2)
  if (nrow(censThoriz)) {
    Surv.censThoriz <- if (is_counting) {
      survfitJM(object, newdata = censThoriz, idVar = idVar, M = M,
                survTimes = Thoriz, last.time = rep(Tstart, length(idcens)),
                LeftTrunc_var = all.vars(TermsT)[1L])
    } else {
      survfitJM(object, newdata = censThoriz, idVar = idVar, M = M,
                survTimes = Thoriz, last.time = rep(Tstart, length(idcens)))
    }
    tt <- Time[indCens]
    weights <- if (is_counting) {
      survfitJM(object, newdata = censThoriz, idVar = idVar, M = M,
                survTimes = Thoriz, last.time = tt[!duplicated(censThoriz[[idVar]])],
                LeftTrunc_var = all.vars(TermsT)[1L])
    } else {
      survfitJM(object, newdata = censThoriz, idVar = idVar, M = M,
                survTimes = Thoriz, last.time = tt[!duplicated(censThoriz[[idVar]])])
    }
    Surv.censThoriz <- sapply(Surv.censThoriz$summaries, "[", 2)
    weights <- sapply(weights$summaries, "[", 2)
  } else {
    Surv.censThoriz <- weights <- NA
  }
  
  # until here I copy paste the prederrJM.mv function from the JMbayes package
  
  
  statusCens <- rbinom(n = length(weights), size = 1, prob = weights)
  matCens <- data.frame(pred_prob = Surv.censThoriz, status = statusCens)
  
  # New part for calibration plot
  matDead <- data.frame(pred_prob = 1 - Surv.deadThoriz, status = 1)
  matAlive <- data.frame(pred_prob = 1 - Surv.aliveThoriz, status = 0)
  if (include.cens == TRUE) {
    mat <- rbind(matDead, matAlive, matCens)
  } else {
    mat <- rbind(matDead, matAlive)
  }
  
  
  # Logistic calibration curve
  y <- mat$status
  p <- mat$pred_prob
  f.recal <- glm(y ~ p, data = mat, family = binomial)
  logit <- seq(-7, 7, length = 200)
  prob <- plogis(logit)
  newdata = list(p = prob)
  log.cal.curve <- predict(f.recal, newdata = newdata, type = "response")
  log.cal.curve <- list(x = prob, y = log.cal.curve)
  # Harrell code - logistic calibration curve
  # y <- mat$status
  # p <- mat$pred_prob
  # logit <- qlogis(p)
  # i <- !is.infinite(logit)
  # f.recal <- lrm.fit(logit[i], y = y[i])
  # logit <- seq(-7, 7, length = 200)
  # prob <- plogis(logit)
  # pred.prob <- f.recal$coef[1] + f.recal$coef[2] * logit 
  # pred.prob2 <- plogis(pred.prob)
  
  # Harrel code - optionally a smooth nonparametric fit
  sm.nonpar.fit <- lowess(mat$pred_prob, mat$status, iter = 0)
  
  # Logistic calibration curve with splines terms
  y <- mat$status
  p <- mat$pred_prob
  f.recal <- glm(y ~ ns(p, 2), data = mat, family = binomial)
  logit <- seq(-7, 7, length = 200)
  prob <- plogis(logit)
  newdata = list(p = prob)
  log.cal.curve.poly <- predict(f.recal, newdata = newdata, type = "response")
  log.cal.curve.poly <- list(x = prob, y = log.cal.curve.poly)
  
  out <- list(log.cal.curve = log.cal.curve,
              log.cal.curve.poly = log.cal.curve.poly,
              sm.nonpar.fit = sm.nonpar.fit)
  
  if (pl == TRUE) {
    plot(0.5, 0.5, xlim = c(0, 1), ylim = c(0, 1), type = "n",
         xlab = "Predicted Probability",
         ylab = "Actual Probability",)
    abline(0, 1, lwd = 6, col = gray(0.85))
    abline(v = max(p), col = 2, lty = 2)
    lines(log.cal.curve, lty = 1, lwd = 2)
    lines(log.cal.curve.poly, lty = 2)
    lines(sm.nonpar.fit, lty = 3)
    leg <- "Ideal"
    leg <- c(leg, "Logistic calibration")
    leg <- c(leg, "Logistic calibration (splines)")
    leg <- c(leg, "Nonparametric")
    leg <- c(leg, "Max predicted prob")
    lp = c(-0.05, 1.05)
    lp = list(x = lp[1], y = lp[2])
    legend(lp, leg, lty = c(1,1,2,3,6), cex = 1, y.intersp = 0.5,
           lwd = c(6,2,1,1,1), col = c(gray(0.85), 1, 1, 1, 2), bty = "n")
  }
  out
  # CHECK plot with val.prob(mat$pred_prob, mat$status)
  # Still to check: polynomial logistic calibration - not sure
}
