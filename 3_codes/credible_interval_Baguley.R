lm.ci <- function(data.frame, conf.level = 0.95, difference = FALSE) {
    #loftus-masson within-subject CIs
    k = ncol(data.frame)
    n <- nrow(data.frame)
    df.stack <- stack(data.frame)
    require(nlme)
    parts <- rep(1:n, k)
    root.ms.error <- lme(values ~ 0 + ind, random = ~1 | parts, cbind(parts, 
                                                                      df.stack))[[6]]
    detach(package:nlme)
    mean.mat <- matrix(NA, k, 1)
    ci.mat <- matrix(NA, k, 2)
    if (difference == TRUE) diff.factor = 2^0.5/2 else diff.factor = 1
    moe <- root.ms.error/n^0.5 * qt(1 - (1 - conf.level)/2, (n - 1) * (k - 
                                                                           1)) * diff.factor
    for (i in 1:k) mean.mat[i,] <- mean(as.matrix(data.frame[,i]))
    for (i in 1:k) {
        ci.mat[i, 1] <- mean.mat[i] - moe
        ci.mat[i, 2] <- mean.mat[i] + moe
    }
    dimnames(ci.mat) <- list(names(data.frame), c("lower", "upper"))
    ci.mat
}

bs.ci <- function(data.frame, conf.level = 0.95, difference = FALSE) {
    # between-subject CIs
    k = ncol(data.frame)
    n <- nrow(data.frame)
    df.stack <- stack(data.frame)
    group.means <- mean(as.matrix(data.frame), na.rm = TRUE)
    if (difference == TRUE) {
        ci.mat <- (confint(lm(values ~ 0 + ind, df.stack)) - group.means) * 
            2^0.5/2 + group.means
    } else {
        ci.mat <- confint(lm(values ~ 0 + ind, df.stack))
    }
    dimnames(ci.mat) <- list(names(data.frame), c("lower", "upper"))
    ci.mat
}

cm.ci.u <- function(data.frame, conf.level = 0.95) {
    #cousineau uncorrected within-subject CIs
    k = ncol(data.frame)
    n <- nrow(data.frame)
    df.stack <- stack(data.frame)
    index <- rep(1:n, k)
    p.means <- tapply(df.stack$values, index, mean)
    norm.df <- data.frame - p.means + (sum(data.frame)/(n * k))
    output <- matrix(NA, k, 2)
    dimnames(output) <- list(names(data.frame), c("lower", "upper"))
    for (i in 1:k) output[i, ] <- t.test(norm.df[i], conf.level = conf.level)$conf.int[1:2]
    output
}

cm.ci <- function(data.frame, conf.level = 0.95, difference = TRUE) {
    #cousineau-morey within-subject CIs
    k = ncol(data.frame)
    if (difference == TRUE) diff.factor = 2^0.5/2 else diff.factor = 1
    n <- nrow(data.frame)
    df.stack <- stack(data.frame)
    index <- rep(1:n, k)
    p.means <- tapply(df.stack$values, index, mean)
    norm.df <- data.frame - p.means + (sum(data.frame)/(n * k))
    t.mat <- matrix(NA, k, 1)
    mean.mat <- matrix(NA, k, 1)
    for (i in 1:k) t.mat[i, ] <- t.test(norm.df[i])$statistic[1]
    for (i in 1:k) mean.mat[i, ] <- mean(norm.df[,i])
    c.factor <- (k/(k - 1))^0.5
    moe.mat <- mean.mat/t.mat * qt(1 - (1 - conf.level)/2, n - 1) * c.factor * 
        diff.factor
    ci.mat <- matrix(NA, k, 2)
    dimnames(ci.mat) <- list(names(data.frame), c("lower", "upper"))
    for (i in 1:k) {
        ci.mat[i, 1] <- mean.mat[i] - moe.mat[i]
        ci.mat[i, 2] <- mean.mat[i] + moe.mat[i]
    }
    ci.mat
}

ml.ci <- function(data.frame, conf.level = 0.95, cov.matrix = "unstructured") {
    # CI based on multilevel model with covariance matrix unstructured or
    #   constrained to compound symmetry
    k = ncol(data.frame)
    n.parts <- nrow(data.frame)
    data.long <- reshape(as.data.frame(data.frame), idvar = "id", direction = "long", varying = 1:k, 
                         v.names = "dv")
    require(nlme)
    if (cov.matrix == "comp.symm") {
        ml.mod <- lme(dv ~ 0 + factor(time), random = ~1 | id, na.action = na.omit, 
                      data.long)
    } else {
        ml.mod <- lme(dv ~ 0 + factor(time), random = ~0 + factor(time) | 
                          id, na.action = na.omit, data.long)
    }
    detach(package:nlme)
    require(gmodels)
    ci.mat <- ci(ml.mod, confidence = conf.level)[, 2:3]
    detach(package:gmodels)
    ci.mat
}

plot.wsci <- function(data.frame, conf.level = 0.95, type = "cm", 
                      difference = TRUE, cov.matrix = "unstructured", xlab = NULL, ylab = NULL, 
                      level.labels = NULL, main = NULL, pch = 23, ylim = c(min.y, max.y),
                      line.width= c(1, 0)) {
    # plot within-subject CIs by various methods
    k = ncol(data.frame)
    if (type == "cm") 
        ci.mat <- cm.ci(data.frame, conf.level, difference = difference)
    if (type == "uncorrected") 
        ci.mat <- cm.ci.u(data.frame, conf.level)
    if (type == "lm") 
        ci.mat <- lm.ci(data.frame, conf.level, difference = difference)
    if (type == "bs") 
        ci.mat <- bs.ci(data.frame, conf.level, difference = difference)
    if (type == "ml") 
        ci.mat <- ml.ci(data.frame, conf.level, cov.matrix = cov.matrix)
    moe.y <- max(ci.mat) - min(ci.mat)
    min.y <- min(ci.mat) - moe.y/3
    max.y <- max(ci.mat) + moe.y/3
    means <- mean(data.frame, na.rm = TRUE)
    if (missing(xlab)) 
        xlab <- "levels"
    if (missing(ylab)) 
        ylab <- "Confidence interval for mean"
    if (missing(level.labels) == FALSE) 
        names(data.frame) <- level.labels
    plot(0, 0, ylim = ylim, xaxt = "n", xlim = c(0.7, k + 0.3), xlab = xlab, 
         ylab = ylab, main = main)
    points(means, pch = pch, bg = "black")
    index <- 1:k
    segments(index, ci.mat[, 1], index, ci.mat[, 2], lwd = line.width[1])
    segments(index - 0.02, ci.mat[, 1], index + 0.02, ci.mat[, 1], lwd = line.width[2])
    segments(index - 0.02, ci.mat[, 2], index + 0.02, ci.mat[, 2], lwd = line.width[2])
    axis(1, 1:k, labels = names(data.frame))
}

two.tiered.ci <- function(data.frame, conf.level = 0.95, cov.matrix = "unstructured", 
                          difference = TRUE, level.labels = NULL, xlab = NULL, ylab = NULL, main = NULL, 
                          pch = 19, pch.cex = 1.4, text.cex = 1.2, ylim = c(min.y, max.y), grid = FALSE,
                          line.width= c(1.5, 1.5)) {
    # plot two tiered CI with ml approach for outer tier and cm approach for inner tier
    k = ncol(data.frame)
    ci.outer <- ml.ci(data.frame, conf.level = conf.level, cov.matrix = cov.matrix)
    moe.y <- max(ci.outer) - min(ci.outer)
    min.y <- min(ci.outer) - moe.y/3
    max.y <- max(ci.outer) + moe.y/3
    means <- mean(data.frame)
    if (missing(xlab)) 
        xlab <- "levels"
    if (missing(ylab)) 
        ylab <- "Confidence interval for mean"
    if (missing(level.labels) == FALSE) 
        names(data.frame) <- level.labels
    plot(0, 0, ylim = ylim, xaxt = "n", xlim = c(0.7, k + 0.3), xlab = xlab, 
         ylab = ylab, main = main, cex.lab = text.cex)
    if (grid == TRUE) 
        grid()
    points(means, pch = pch, bg = "black", cex = pch.cex)
    index <- 1:k
    segments(index, ci.outer[, 1], index, ci.outer[, 2], lwd = line.width[1])
    axis(1, index, labels = names(data.frame))
    row.names(ci.outer) <- names(data.frame)
    if (cov.matrix == "comp.symm") 
        ci.inner <- lm.ci(data.frame, conf.level, difference = difference)
    else ci.inner <- cm.ci(data.frame, conf.level, difference = difference)
    segments(index - 0.02, ci.inner[, 1], index + 0.02, ci.inner[, 1], lwd = line.width[1])
    segments(index - 0.02, ci.inner[, 2], index + 0.02, ci.inner[, 2], lwd = line.width[1])
}

cm.ci.mixed <- function(data.frame, group.var = "last", conf.level = 0.95, 
                        difference = TRUE) {
    #cousineau-morey within-subject CIs for mixed design
    k = ncol(data.frame)
    if (difference == TRUE) diff.factor = 2^0.5/2 else diff.factor = 1
    if (group.var == "last") {
        within <- 1:(k - 1)
        data.frame[k] <- unclass(data.frame[[k]])[1:nrow(data.frame)]
    } else {
        within <- 2:k
        data.frame[1] <- unclass(data.frame[[1]])[1:nrow(data.frame)]
    }
    if (group.var == "last") n.groups <- nlevels(factor(data.frame[[k]])) else n.groups <- nlevels(factor(data.frame[[1]]))
    c.factor <- ((k - 1)/(k - 2))^0.5
    ci.list <- vector("list", n.groups)
    for (i in 1:n.groups) {
        if (group.var == "last") data <- subset(data.frame, data.frame[k] == i)[within] else data <- subset(data.frame, data.frame[1] == i)[within]
        p.means <- mean(as.matrix(t(data)))
        norm.dat <- data - p.means + (mean(p.means))
        t.mat <- matrix(NA, k - 1, 1)
        mean.mat <- matrix(NA, k - 1, 1)
        ci.mat <- matrix(NA, k - 1, 2)
        for (j in 1:(k - 1)) t.mat[j, ] <- t.test(norm.dat[j])$statistic[1]
        for (j in 1:(k - 1)) mean.mat[j, ] <- mean(norm.dat[,j])
        n <- nrow(data)
        moe.mat <- mean.mat/t.mat * qt(1 - (1 - conf.level)/2, n - 1) * c.factor * 
            diff.factor
        for (j in 1:(k - 1)) {
            ci.mat[j, 1] <- mean.mat[j] - moe.mat[j]
            ci.mat[j, 2] <- mean.mat[j] + moe.mat[j]
        }
        dimnames(ci.mat) <- list(names(data), c("lower", "upper"))
        ci.list[[i]] <- ci.mat
    }
    ci.list
}

ml.ci.mixed <- function(data.frame, group.var = "last", conf.level = 0.95, 
                        cov.matrix = "unstructured") {
    # mixed CI based on multilevel model unstructured or constrained to
    #   compound symmetry
    k = ncol(data.frame)
    n.parts <- nrow(data.frame)
    if (group.var == "last") 
        within <- 1:(k - 1)
    else within <- 2:k
    data.long <- reshape(data.frame, idvar = "id", direction = "long", varying = within, 
                         v.names = "dv")
    group <- factor(data.long[[1]])
    n.groups <- nlevels(group)
    require(nlme)
    if (cov.matrix == "within.group.cs") 
        ml.mod <- lme(dv ~ 0 + factor(time):factor(group), random = ~0 + 
                          factor(time) | id, na.action = na.omit, cbind(group, data.long))
    if (cov.matrix == "unstructured") 
        ml.mod <- lme(dv ~ 0 + factor(time):factor(group), random = list(id = pdDiag(form = ~0 + 
                                                                                         factor(group)), id = ~0 + factor(time)), na.action = na.omit, 
                      cbind(group, data.long))
    if (cov.matrix == "comp.symm") 
        ml.mod <- lme(dv ~ 0 + factor(time):factor(group), random = ~1 | 
                          id, na.action = na.omit, cbind(group, data.long))
    require(gmodels)
    detach(package:nlme)
    ci.mat <- ci(ml.mod, confidence = conf.level)[, 2:3]
    group.id <- rep(1:n.groups, rep(k - 1, n.groups))
    output <- cbind(group.id, ci.mat)
    detach(package:gmodels)
    output
}

two.tiered.mixed <- function(data.frame, group.var = "last", conf.level = 0.95, 
                             cov.matrix = "unstructured", difference = TRUE, lines = FALSE, level.labels = NULL, 
                             xlab = NULL, ylab = NULL, main = NULL, pch = c(21:25, 1:3), pch.cex = 1.2, 
                             pch.col = c(3:9), text.cex = 1.2, ylim = c(min.y, max.y), grid = FALSE, jitter = NULL,
                             group.labels = c(1:8), leg.loc = c(1, min(ci.outer)), line.width= c(1.25, 1.5, 1)) {
    # plot two tiered mixed CI with ml approach for outer tier and cm
    #   approach for inner tier
    k = ncol(data.frame)
    if (group.var == "last") 
        n.groups <- nlevels(factor(data.frame[[k]]))
    else n.groups <- nlevels(factor(data.frame[[1]]))
    ci.outer <- ml.ci.mixed(data.frame, group.var = group.var, conf.level = conf.level, 
                            cov.matrix = cov.matrix)
    ci.inner <- cm.ci.mixed(data.frame, group.var = group.var, conf.level, 
                            difference = difference)
    ci.group.outer <- matrix(NA, k - 1, 2)
    ci.group.inner <- matrix(NA, k - 1, 2)
    index <- 1:(k - 1)
    moe.y <- max(ci.outer) - min(ci.outer)
    min.y <- min(ci.outer) - moe.y/3
    max.y <- max(ci.outer) + moe.y/3
    if (missing(xlab)) 
        xlab <- "levels"
    if (missing(ylab)) 
        ylab <- "Confidence interval for mean"
    if (missing(level.labels) == FALSE) 
        names(data.frame) <- level.labels
    plot(0, 0, ylim = ylim, xaxt = "n", xlim = c(0.7, k - 1 + 0.3), xlab = xlab, 
         ylab = ylab, main = main, cex.lab = text.cex)
    if (grid == TRUE) 
        grid()
    axis(1, index, labels = rownames(ci.inner[[1]]))
    if (missing(jitter)) 
        jitter <- scale(1:n.groups, scale = FALSE)/(n.groups * 1.5)
    for (i in 1:n.groups) {
        ci.group.outer <- subset(as.data.frame(ci.outer), group.id == i)[2:3]
        ci.group.inner <- ci.inner[[i]]
        means <- (ci.group.inner[, 1] + ci.group.inner[, 2])/2
        points(index + jitter[i], means, pch = pch[i], bg = pch.col[i], cex = pch.cex)
        segments(index + jitter[i], ci.group.outer[, 1], index + jitter[i], 
                 ci.group.outer[, 2], lwd = line.width[1])
        segments(index - 0.02 + jitter[i], ci.group.inner[, 1], index + 0.02 + 
                     jitter[i], ci.group.inner[, 1], lwd = line.width[2])
        segments(index - 0.02 + jitter[i], ci.group.inner[, 2], index + 0.02 + 
                     jitter[i], ci.group.inner[, 2], lwd = line.width[2])
        if (lines == TRUE) 
            lines(index + jitter[i], means, lty = i + 1, lwd = line.width[3])
    }
    legend(leg.loc[1], leg.loc[2], legend = group.labels[1:n.groups], pch = pch[1:n.groups], 
           lty = 2:(n.groups + 1), pt.bg = pch.col[1:n.groups], bty = "n", horiz = TRUE, lwd = line.width[3])
} 