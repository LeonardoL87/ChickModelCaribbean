library(grDevices)
spectra <- read.csv("spectra.csv")
#head(spectra)
spectra$Frequency <- as.numeric(sub(",", ".",  sub(".", "", spectra$Frequency, fixed = TRUE), fixed = TRUE))
spectra$dB.Spec. <- as.numeric(sub(",", ".",  sub(".", "", spectra$dB.Spec., fixed = TRUE), fixed = TRUE))
spectra$Frequency.1 <- as.numeric(sub(",", ".",  sub(".", "", spectra$Frequency.1, fixed = TRUE), fixed = TRUE))
spectra$dB.Spec..1 <- as.numeric(sub(",", ".",  sub(".", "", spectra$dB.Spec..1, fixed = TRUE), fixed = TRUE))
pdf('spectra_obs_vs_fitted.pdf', height = 7 , width = 5, bg = "white")
plot(spectra$Frequency, spectra$dB.Spec., col="blue", type="l", lwd=3,
     xlab = "Frequency", ylab = "Amplitud (dB)", main = "Spectra")
lines(spectra$Frequency.1, spectra$dB.Spec..1, col="red", type="l", lwd=3,
     xlab = "Frequency", ylab = "Amplitud (dB)")
#pdf(filename = "spectra_obs_vs_fitted.pdf", height = 295, width = 300, bg = "white")
#dev.off()