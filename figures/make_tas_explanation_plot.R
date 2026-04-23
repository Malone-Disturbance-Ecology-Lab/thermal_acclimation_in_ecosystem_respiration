tot <- read.csv("data/outcome_temp.csv")
dir <- read.csv("data/outcome_temp_water_gpp.csv")

m <- merge(
  tot[, c("site_ID", "TAS", "TASp")],
  dir[, c("site_ID", "TAS", "TASp")],
  by = "site_ID",
  suffixes = c("_tot", "_direct")
)
m$TAS_app <- m$TAS_tot - m$TAS_direct

siteyear_tot <- read.csv("data/outcome_siteyear_temp.csv")
siteyear_dir <- read.csv("data/outcome_siteyear_temp_water_gpp.csv")
m_siteyear <- merge(
  siteyear_tot[, c("site_ID", "growing_year", "window", "lnRatio")],
  siteyear_dir[, c("site_ID", "growing_year", "window", "lnRatio")],
  by = c("site_ID", "growing_year", "window"),
  suffixes = c("_tot", "_dir")
)

bysite_corr <- do.call(
  rbind,
  lapply(split(m_siteyear, m_siteyear$site_ID), function(g) {
    if (nrow(g) < 5) {
      return(NULL)
    }
    data.frame(
      site_ID = g$site_ID[1],
      corr_lnRatio = cor(g$lnRatio_tot, g$lnRatio_dir, use = "complete.obs")
    )
  })
)

cor_td <- cor(m$TAS_tot, m$TAS_direct, use = "complete.obs")
cor_ta <- cor(m$TAS_tot, m$TAS_app, use = "complete.obs")
var_tot <- var(m$TAS_tot, na.rm = TRUE)
var_dir <- var(m$TAS_direct, na.rm = TRUE)
var_app <- var(m$TAS_app, na.rm = TRUE)
cov_ad <- cov(m$TAS_app, m$TAS_direct, use = "complete.obs")
share_app <- (var_app + cov_ad) / var_tot
share_dir <- (var_dir + cov_ad) / var_tot
rmse_td <- sqrt(mean((m$TAS_direct - m$TAS_tot)^2, na.rm = TRUE))
same_sign <- mean(m$TAS_tot * m$TAS_direct > 0, na.rm = TRUE)
median_within_corr <- median(bysite_corr$corr_lnRatio, na.rm = TRUE)

png(
  filename = "figures/tas_explanation.png",
  width = 1800,
  height = 1600,
  res = 180
)

op <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), mar = c(4.8, 4.8, 4.2, 1.2), oma = c(0, 0, 2.2, 0), cex.main = 0.9)

# Panel A: estimator logic
plot.new()
plot.window(xlim = c(0, 10), ylim = c(0, 10))
rect(0.8, 7.2, 4.3, 9.15, col = "#dceefb", border = "#3b6ea8", lwd = 2)
text(2.55, 8.6, "Total estimator", cex = 1.1, font = 2)
text(2.55, 8.0, "Fit ER(TS) in each year-window", cex = 0.8)
text(2.55, 7.55, "Then regress: lnRatio ~ TS + window", cex = 0.8)

rect(5.7, 7.2, 9.2, 9.15, col = "#e5f5e0", border = "#4a8f4d", lwd = 2)
text(7.55, 8.65, "Direct estimator", cex = 1.1, font = 2)
text(7.45, 8.0, "Fit ER(TS, SWC, daytime NEE)", cex = 0.8)
text(7.45, 7.55, "Evaluate at fixed window references", cex = 0.8)

arrows(2.45, 7.35, 2.45, 5.8, length = 0.08, lwd = 2, col = "#3b6ea8")
arrows(7.55, 7.35, 7.55, 5.8, length = 0.08, lwd = 2, col = "#4a8f4d")

rect(1.1, 4.05, 4.0, 5.7, col = "#dceefb", border = "#3b6ea8", lwd = 2)
text(2.45, 5.1, "Slope = TAS_tot", cex = 1.05, font = 2)
text(2.45, 4.45, "Reduced-form TS slope", cex = 0.9)

rect(6.0, 4.05, 8.9, 5.7, col = "#e5f5e0", border = "#4a8f4d", lwd = 2)
text(7.55, 5.1, "Slope = TAS_direct", cex = 1.05, font = 2)
text(7.55, 4.45, "Conditional TS slope", cex = 0.9)

arrows(3.9, 4.8, 6.05, 4.8, length = 0.08, lwd = 2, col = "#555555")
text(4.98, 5.25, "difference", cex = 0.9)

rect(2.8, 1.0, 7.2, 2.8, col = "#fff3cd", border = "#c89f2d", lwd = 2)
text(5.0, 2.2, "TAS_app = TAS_tot - TAS_direct", cex = 1.15, font = 2)
text(5.0, 1.55, "Residual contrast, not a separately identified mediation effect", cex = 0.8)
title("A. What the code estimates", line = 0.8, font.main = 2)

# Panel B: total vs direct
xlim_td <- range(c(m$TAS_tot, m$TAS_direct), na.rm = TRUE)
plot(
  m$TAS_tot, m$TAS_direct,
  pch = 19, cex = 0.7, col = rgb(44, 127, 184, 140, maxColorValue = 255),
  xlab = "TAS_tot", ylab = "TAS_direct",
  xlim = xlim_td, ylim = xlim_td
)
abline(0, 1, lty = 2, lwd = 2, col = "#7f7f7f")
abline(lm(TAS_direct ~ TAS_tot, data = m), col = "#08519c", lwd = 2)
usr <- par("usr")
text(
  usr[1] + 0.06 * diff(usr[1:2]),
  usr[4] - 0.04 * diff(usr[3:4]),
  labels = paste(
    sprintf("cor = %.2f, RMSE = %.3f", cor_td, rmse_td),
    sprintf("same sign = %.0f%%", 100 * same_sign),
    sprintf("var(TAS_tot) = %.4f", var_tot),
    sprintf("var(TAS_direct) = %.4f", var_dir),
    sep = "\n"
  ),
  adj = c(0, 1),
  cex = 0.82
)
title("B. Direct tracks total, but is not just a noisier copy", line = 0.8, font.main = 2)

# Panel C: total vs apparent
plot(
  m$TAS_tot, m$TAS_app,
  pch = 19, cex = 0.7, col = rgb(217, 95, 2, 140, maxColorValue = 255),
  xlab = "TAS_tot", ylab = "TAS_app = TAS_tot - TAS_direct"
)
abline(h = 0, v = 0, lty = 3, col = "#7f7f7f")
abline(lm(TAS_app ~ TAS_tot, data = m), col = "#b35806", lwd = 2)
usr <- par("usr")
text(
  usr[1] + 0.06 * diff(usr[1:2]),
  usr[4] - 0.04 * diff(usr[3:4]),
  labels = paste(
    sprintf("cor = %.2f", cor_ta),
    sprintf("var(TAS_app) = %.4f", var_app),
    "Residual still has substantial spread",
    "So it is not negligible or purely noise",
    sep = "\n"
  ),
  adj = c(0, 1),
  cex = 0.82
)
title("C. The apparent component is substantial", line = 0.8, font.main = 2)

# Panel D: shares + within-site correlation
bp <- barplot(
  c(share_dir, share_app),
  names.arg = c("Direct share", "Apparent share"),
  ylim = c(0, 1.28),
  col = c("#66c2a4", "#fdb462"),
  border = NA,
  ylab = "Share of TAS_tot variance\n(var + covariance term) / var(TAS_tot)"
)
text(bp, c(share_dir, share_app) + 0.05, labels = sprintf("%.2f", c(share_dir, share_app)), cex = 1)
abline(h = 1, lty = 2, col = "#7f7f7f")
usr <- par("usr")
text(mean(usr[1:2]), 1.16,
     labels = sprintf("Median within-site corr of lnRatio_tot vs lnRatio_direct = %.2f", median_within_corr),
     cex = 0.8)
text(mean(usr[1:2]), 1.09,
     labels = "Window fixed effects absorb simple between-window level differences.",
     cex = 0.76)
text(mean(usr[1:2]), 1.02,
     labels = "What remains is a residual contrast after conditioning,",
     cex = 0.76)
text(mean(usr[1:2]), 0.96,
     labels = "not clean mediation identification.",
     cex = 0.76)
title("D. Why the residual should be described carefully", line = 0.8, font.main = 2)

mtext("Thermal response decomposition in the current workflow", outer = TRUE, cex = 1.35, font = 2)
par(op)
dev.off()
