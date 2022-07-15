######################################################################
#This programme calculates####
#the structure function given the RA and Dec and P% and position######
#angle. 
######################################################################
import math
import numpy as np
import matplotlib.pyplot as plt
# from strcompress import *
from scipy.optimize import curve_fit
import matplotlib.ticker as tik
import pylab as p
######################################################################
############Calling the columns of the file in Python#################
######################################################################
#n1, ra1, dec1, P, eP, PA, ePA = np.loadtxt('Pol_ra_dec_vectors_rest_1.dat', unpack=True)
#num_lines = sum(1 for line in open('Pol_ra_dec_vectors_rest_1.dat'))

#n1, ra1, dec1, P, eP, PA, ePA = np.loadtxt('Pol_ra_dec_vectors_UptoFields5_1.dat', unpack=True)
#num_lines = sum(1 for line in open('Pol_ra_dec_vectors_UptoFields5_1.dat'))

n1, ra1, dec1, P, eP, PA, ePA = np.loadtxt('L1415_pol_data.dat', unpack=True)
num_lines = sum(1 for line in open('L1415_pol_data.dat'))


print('Total number of entries=',num_lines)
print('Total number of calculated entries=',num_lines*(num_lines-1)/2)
####################Making the calculations###########################
count = 0
f1 = open("diff_l1415.dat","w")
while (count <= num_lines-1):
  for i in range(0, num_lines):
    if count < i:
      diff_PA	= math.pow((PA[count] - PA[i]),2)
      diff_ePA  = (ePA[count]**2 + ePA[i]**2) 
      ra_diff	= ra1[count] - ra1[i]
      dec_diff	= dec1[count] - dec1[i]
      sep_l	= math.sqrt(math.pow(ra_diff,2) + math.pow(dec_diff,2))
      sep_li_pc = 250*3600*5*math.pow(10,-6)*sep_l
      f1.write('{0:d} \t {1:d} \t {2:f} \t {3:f} \t {4:f} \t {5:f} \t {6:f} \t {7:f} \n'.format(count, i, diff_PA, diff_ePA, sep_l, sep_li_pc, ra1[i], dec1[i]))
  count += 1
f1.close()
print ('Total number of calculated entries=',num_lines*(num_lines-1)/2)
##################Called the file --> newfile#########################
count1 ,i1 ,diff_PA1, diff_ePA1, sep_l1, sep_li_pc1, ra2, dec2 = np.loadtxt('diff_l1415.dat', unpack=True)

#plt.hist(sep_li_pc1, bins = 100)
#plt.show()

for l in range(len(sep_li_pc1)):
  if sep_li_pc1[l] <= 0.01:
    print(sep_li_pc1[l], diff_PA1[l])

min_val  = 0.05#min(sep_li_pc1)
max_val  = np.max(sep_li_pc1)
#min_val  = 0.1#min(sep_l1)
#max_val  = max(sep_l1)

ratio = np.log(max_val/min_val)
x = 1.25

bins = np.int(ratio/np.log(x))
print (bins)
'''
ratio = (max_val/min_val)
x = 1.2

bins = np.int(ratio/(x))
print (bins)
'''
f1 = open("test.dat","w")
for i in range(bins):
  range_1 = min_val * math.pow(x,i)
  f1.write('{0:f} \n'.format(range_1))
  print(range_1)
f1.close()
range_= np.loadtxt('test.dat', unpack=True)
num_lines2 = sum(1 for line in range_)
print (min_val, max_val)
######################################################################
f2 = open("mean_PositionA","w")
count1 = 0
while (count1 < num_lines2 - 1):
  filename = ('temp'+str(count1))
  f3 = open(filename,"w")
  for j1 in range(len(sep_li_pc1)):
    if sep_li_pc1[j1] >= range_[count1] and sep_li_pc1[j1] < range_[count1+1]:
      mid_range = (range_[count1] + range_[count1+1])/2
      f3.write('{0:f} \t {1:f} \t {2:f} \t {3:f} \t {4:f} \t {5:f} \n'.format(diff_PA1[j1], diff_ePA1[j1], sep_li_pc1[j1], sep_l1[j1], ra2[j1], dec2[j1]))
  f3.close()
  diff_PA2, diff_ePA2 = np.loadtxt(filename, unpack=True, usecols=[0,1])
  num_lines3 = len(diff_PA2)
  #sum_PA = sum(diff_PA2)/num_lines3
  #sum_ePA = sum(diff_ePA2)/num_lines3
  #std_PA = sqrt(sum_PA - sum_ePA)
  std_PA = (np.mean(diff_PA2)) - (np.mean(diff_ePA2))
  mean_PA = np.sqrt(std_PA)
  #std_PA = (np.mean(diff_PA2))
  #mean_PA = np.sqrt(std_PA)
  #err_PA1 = (math.sqrt(sum((mean_PA-np.sqrt(diff_PA2))**2))/num_lines3)
  #err_PA = err_PA1/math.sqrt(num_lines3)
  a = np.sqrt(diff_ePA2)
  err_PA1 = np.std(a)
  err_PA = err_PA1/math.sqrt(num_lines3)

  #mean_PA = math.sqrt(mean(diff_PA2))
  print (num_lines3, err_PA1, err_PA, mean_PA)
  f2.write('{0:f} \t {1:f} \t {2:f} \t {3:f}\n'.format(mean_PA, mid_range, err_PA1, err_PA))
  count1 += 1
f2.close()
######################################################################
mean_PA3, sep_l3, diff_ePA3 = np.loadtxt('mean_PositionA', unpack=True, usecols= [0, 1, 3])
mean_PA31, sep_l31 = np.loadtxt('mean_PositionA1', unpack=True, usecols= [0, 1])
##################Plot################################################

def fitfun(l, m, b):
  y = np.sqrt(b**2 + (m**2) * (l**2))
  return y

popt, pcov = curve_fit(fitfun, sep_l31, mean_PA31)

ymodel = fitfun(sep_l31, popt[0], popt[1])

#ymodel1 = [2.465,2.466,3.0,3.4,4.2,5.0,6.0,7.0,7.5,7.8,7.9,8.0,8.2,8.5,8.6]

b = popt[1]
b1 = popt[1]*math.pi/180

perr = np.sqrt(np.diag(pcov))

Bmag = b1/(math.sqrt(2-(b1**2)))

print(b, b1, Bmag,perr)
######################################################################

ax1 = plt.gca()

#fig, ax = plt.plot(ncols = 1, figsize=(6,6), sharey=True)

ax1.errorbar(sep_l3, mean_PA3, yerr=diff_ePA3, linestyle='', marker='o')
ax1.set_xlim([0.1,1.5])
ax1.set_xscale('log')
ax1.set_ylim([10,20])
ax1.step(sep_l3, mean_PA3, where = 'mid')
ax1.plot(sep_l31, ymodel, '--')
ax1.set_xlabel('Distance(pc)')
ax1.set_ylabel('Angular Dispersion<($\Delta$ $\phi$)$^{2}$>$^{1/2}$ ($\degree$)')

ax1.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax1.get_xaxis().set_major_formatter(tik.ScalarFormatter())

ax1.text(1, 19, 'L1415')

plt.show()

