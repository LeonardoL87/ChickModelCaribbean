
import numpy as np
import pylab as pl
import pandas as pd
import seaborn

spectra = pd.read_csv("./spectra_corrected.csv")
freq_sim = spectra["Frequency"]
ampli_sim = spectra["dB.Spec."]
freq_obs = spectra["Frequency.1"]
ampli_obs = spectra["dB.Spec..1"]

# plotting figs:
# some configuration first:
xlabelsize = 15
ylabelsize = 15
pl.rcParams['xtick.labelsize'] = xlabelsize
pl.rcParams['ytick.labelsize'] = ylabelsize
hfont = {'fontname':'Helvetica'}
sefont = {'fontname':'Serif'}
pl.figure(figsize=(10,5))

pl.plot(freq_sim,ampli_sim)
