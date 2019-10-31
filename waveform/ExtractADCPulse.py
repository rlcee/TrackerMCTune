import numpy as np

data = np.genfromtxt("FitADCPulse2.dat")
print "LBL ADC pole  : %5.2f +- %5.2f" % (np.mean(data[:,6]),np.std(data[:,6]))

data = np.genfromtxt("FitADCPulse.dat")
print "Fermi ADC pole: %5.2f +- %5.2f" % (np.mean(data[:,6]),np.std(data[:,6]))
