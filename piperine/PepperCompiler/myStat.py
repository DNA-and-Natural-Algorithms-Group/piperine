"""Some useful statistics methods."""
from __future__ import division

import math
from math import sqrt, exp, log
import copy

NaN = float("nan")

## Sample statistics

def mean(xs):
  """Sample mean."""
  if len(xs) > 0:
    return sum(xs) / len(xs)
  else:
    return NaN

def var(xs):
  """Sample variance."""
  if len(xs) > 1:
    m = mean(xs)
    xs2 = [(x-m)**2 for x in xs]
    return sum(xs2) / (len(xs) - 1)
  else:
    return NaN

def stddev(xs):
  """Sample standard deviation."""
  return sqrt(var(xs))

def moment(k, xs):
  """Sample n-th central moment. Note: not unbiased."""
  if len(xs) > 0:
    m = mean(xs)
    xs_k = [(x-m)**k for x in xs]
    return sum(xs_k) / len(xs)
  else:
    return NaN

def skewness(xs):
  """Sample skewness."""
  return moment(3, xs) / var(xs)**(3/2)

def kurtosis(xs):
  """Sample kurtosis."""
  return moment(4, xs) / var(xs)**2 - 3

def _get_index(xs, n):
  """Get a fractionally indexed item by interpolation."""
  if len(xs) == 0:
    return NaN
  assert 0 <= n <= len(xs) - 1
  if n == len(xs) - 1:
    return xs[-1]
  i = int(n // 1) # Integral part
  f = n % 1  # Fractional part
  
  return (1-f) * xs[i] + f * xs[i+1]

def quantiles(xs, *quantiles):
  """Find a set of quantiles in data."""
  # NOTE: Will be inefficient for large len(xs)
  xs = copy.copy(xs)
  xs.sort()
  
  return [_get_index(xs, quant * (len(xs) - 1)) for quant in quantiles]

def median(xs):
  """Sample median."""
  # NOTE: Will be inefficient for large len(xs)
  return quantiles(xs, 0.5)

def quartiles(xs):
  """Divisions between quartiles."""
  return quantiles(xs, 0.25, 0.5, 0.75)


## Statistical Distributions

# Parameters in Lin approximation
a1 = -0.9911; b1 =  0.8055
a2 = -0.6763; b2 = -1.2451

def chi_dist(deg, x):
  """Returns the cumulative distribution of chi squared up to x.
  Currently, it is an approximation from Lin-1988 <http://www.jstor.org/stable/2348373>"""
  
  z = sqrt(x) - sqrt(deg)
  
  if z <= 0:
    return 1 - exp(b1 * z + a1 * z**2) / 2
  else:
    return exp(b2 * z + a2 * z**2) / 2

def chi_inv(deg, p):
  """The inverse cumulative distribution of chi squared 
  for 'deg' degrees of freedom and 'p' cumulative probability.
  Currently, it is an approximation from Lin-1988 <http://www.jstor.org/stable/2348373>"""
  
  if p >= 0.5:
    c = -log(2 * (1-p))
    z = (-b1 + sqrt(b1**2 - 4*a1*c)) / (2*a1)
  else:
    c = -log(2 * p)
    z = (-b2 - sqrt(b2**2 - 4*a2*c)) / (2*a2)
  
  return (z + sqrt(deg))**2
    

## Exponential Distribution analysis
#
# f(x;m) = 1/m e^{-x/m}
#
# if X ~ f(m) then E[X] = m
# Likewise, the maximum likelihood estimator for m is mean(X's)

def exp_mean_interval(xs, alpha):
  """
  Returns the 100(1-alpha)% confidence interval for the parameter m (the mean) based on data.
  Based on <http://en.wikipedia.org/wiki/Exponential_distribution#Maximum_likelihood>
  """
  m_hat = mean(xs) # Maximum likelihood estimator for m
  n = len(xs)
  
  chi_low = chi_inv(2*n, alpha/2) # Lower chi squared parameter
  lower_bound = m_hat * 2*n / chi_low
  
  chi_up  = chi_inv(2*n, 1 - alpha/2) # Upper parameter
  upper_bound = m_hat * 2*n / chi_up
  
  return lower_bound, m_hat, upper_bound


## Gamma Distribution analysis

def gamma_mle(xs):
  """
  Maximum likelihood estimation for k, theta for a gamma distribution.
  Based on <http://en.wikipedia.org/wiki/Gamma_distribution#Maximum_likelihood_estimation>
  """
  N = len(xs)
  if N < 1:
    return NaN, NaN
  
  s = math.log(sum(xs)/N) - sum(map(math.log, xs))/N
  if s != 0:
    k_hat = (3 - s + math.sqrt( (s-3)**2 + 24*s )) / (12*s)  # Approx mle for k
  else:
    k_hat = NaN
  
  if k_hat != 0:
    theta_hat = sum(xs) / (N * k_hat)
  else:
    theta_hat = NaN
  
  return k_hat, theta_hat
