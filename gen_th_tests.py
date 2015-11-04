# Just building a script for makin' toeholds

from __future__ import division

import numpy as np
import matplotlib
# matplotlib.use('QT4Agg') 
import matplotlib.pyplot as plt
import re
import logging

import stickydesign as sd
import energyfuncs_james as energetics
from tolds_utils import e_barplot

# The purpose of the following test is to produce a landscape depicting the ave
# rage number of toeholds produced by stickydesign as a function of the maximum 
# energy deviation and maximum spurious interactions
n = 5 
spurmx_vec = np.linspace(1.0,5.0, n)
deviat_vec = np.linspace(0.2,2.0, n)
result_mat = np.zeros((n, n))

ef = energetics.energyfuncs()
ef.targetdG = 7.7
reps = 3
th_e = 7.7 
max_ths = 20
th_len = 7

for i in np.arange(n):
  spurmx = spurmx_vec[i]
  for j in np.arange(n):
    deviat = deviat_vec[j]
    convenience_vec = np.zeros(reps)
    for x in np.arange(reps):
      # Generate three first toeholds
      try:
          ends = sd.easyends('TD', th_len, max_ths, interaction=th_e,
                             fdev=deviat/th_e, alphabet='h', adjs=['h','d'],
                             maxspurious=spurmx/th_e, energyfuncs=ef)
          convenience_vec[x] = len(ends)
      except ValueError as e:
          convenience_vec[x] = 0
    result_mat[i, j] = convenience_vec.mean()
    np.savetxt('available_ths.csv', result_mat, delimiter=',') 
