from __future__ import division

import numpy as np
import re

nt = { 'a': 0, 'c': 1, 'g': 2, 't': 3 }
NT = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 }
NT_to_num = lambda s: np.array([ NT[let] for let in s ])
tops = lambda s: 4*s[:,:-1]+s[:,1:]
comp = lambda s: 3 - s[::-1]

def e_barplot(toe_holds, ef_j, labels = None):
  '''
  This function accepts an endarray object and presents a barplot depicting the 
  binding energies of the sequences in the endarray, returning the values as a v
  ector.'''
  
  import matplotlib
  #matplotlib.use('QT4Agg') 
  import matplotlib.pyplot as plt
  
  if labels is None:
    labels = [spec + ' ' + pos \
                for pos in ['A', 'B', 'C'] \
                  for spec in ['First', 'Second']]
  # Grab energy values
  e_vec = ef_j.matching_uniform(toe_holds)
  e_vec_ext = ef_j.th_external_dG(toe_holds)
  e_vec_int = ef_j.th_internal_dG(toe_holds)
  e_vec_avg = (e_vec_ext + e_vec_int) / 2
  # Plot parameters
  n_seq = np.size(toe_holds, 0)
  inds = np.arange(n_seq) + 1
  wid = 0.35
  # Plot
  fig = plt.figure()
  ax = fig.add_subplot(111)
  bars_ext = ax.bar(inds, e_vec_ext, wid, color='b')
  bars_int = ax.bar(inds+wid, e_vec_int, wid, color='g')
  line_targetdG = ax.plot([0, inds[-1] + 2.5 * wid], \
                          2*[ef_j.targetdG], \
                          color='k')
  # Plot labels
  ax.set_ylabel(r'$\delta G_{37}$')
  ax.set_title('Binding energies of generated toeholds')
  ax.legend( (bars_ext[0], bars_int[0], line_targetdG[0]), \
             ('Internal Context', 'External Context', 'Target dG'),
             loc='lower right')
  ax.set_xticks(inds)
  plt.ylim(0, 9)
  plt.xlim(1 - wid/2, n_seq + 2.5*wid)
  x_ticks = ax.set_xticklabels(labels)
  plt.setp(x_ticks, rotation=45, fontsize=10)

def get_th_from_domains(domains_fn):
  from Test_Designer_modules import Read_Template
  domains = Read_Template(domains_fn)
  f_A = NT_to_num(domains['f_A'])
  f_B = NT_to_num(domains['f_B'])
  f_C = NT_to_num(domains['f_C'])
  s_A = comp(NT_to_num(domains['s_A_comp']))
  s_B = comp(NT_to_num(domains['s_B_comp']))
  s_C = comp(NT_to_num(domains['s_C_comp']))
  first_ths = np.array( (f_A, f_B, f_C))
  secon_ths = np.array( (s_A, s_B, s_C))
  c_col = np.ones((3,1))
  secon_th_col = np.array(secon_ths[:,0], ndmin=2).transpose()
  first_th_col = np.array(first_ths[:,-1], ndmin=2).transpose()
  first_ths = np.concatenate((c_col, first_ths, secon_th_col), 1)
  secon_ths = np.concatenate( (first_th_col, secon_ths, c_col), 1)
  return np.concatenate( (first_ths, secon_ths), 0)

def gen_random_fixed(th_l=7, bm_l=15, cl=2, fname='ran.fixed'):
  import random
  H = ['A', 'T', 'C']
  W = ['A', 'T']
  N = ['A', 'C', 'G', 'T']
  if False:
    H = N
    W = N
  dom_dict = {
              'f_A':[ random.choice(H) for x in range(th_l) ],
              'f_B':[ random.choice(H) for x in range(th_l) ],
              'f_C':[ random.choice(H) for x in range(th_l) ],
              's_A':[ random.choice(H) for x in range(th_l) ],
              's_B':[ random.choice(H) for x in range(th_l) ],
              's_C':[ random.choice(H) for x in range(th_l) ],
              'm_A':[ random.choice(H) for x in range(bm_l - 2) ],
              'm_B':[ random.choice(H) for x in range(bm_l - 2) ],
              'm_C':[ random.choice(H) for x in range(bm_l - 2) ],
              'h_A_p':['C'] + [ random.choice(H) for x in range(bm_l - 3) ],
              'h_A_q':['C'] + [ random.choice(H) for x in range(bm_l - 3) ],
              'h_B_r':['C'] + [ random.choice(H) for x in range(bm_l - 3) ],
              'h_B_s':['C'] + [ random.choice(H) for x in range(bm_l - 3) ],
              'h_C_j':['C'] + [ random.choice(H) for x in range(bm_l - 3) ],
              'h_C_k':['C'] + [ random.choice(H) for x in range(bm_l - 3) ],
              'm_A_cl':[ random.choice(W), 'C' ],
              'm_B_cl':[ random.choice(W), 'C' ],
              'm_C_cl':[ random.choice(W), 'C' ],
              'h_A_p_cl':[ 'C', random.choice(W)],
              'h_A_q_cl':[ 'C', random.choice(W)],
              'h_B_r_cl':[ 'C', random.choice(W)],
              'h_B_s_cl':[ 'C', random.choice(W)],
              'h_C_j_cl':[ 'C', random.choice(W)],
              'h_C_k_cl':[ 'C', random.choice(W)]}
  f = open(fname, 'w')
  for out_dom in dom_trans.keys():
    f.write('sequence ' + out_dom + ' = ' + 
            str.join('',dom_dict[dom_trans[out_dom]]) + '\n')
  f.close()

def verboten_finder(st_file='Oscillator_james_35.st'):
  seq_str = open(st_file, 'r').read()
  verb_strs = 'GGG|GGGG|CCC|CCCC|TTTT|AAAA|GCGC|GGCC|CCGG|CGCG'
  verb_pats = '|[AT]{6}|[GC]{6}|[TC]{6}|[AG]{6}|[TC][AG][TC][AG][TC][AG]|[AG][TC][AG][TC][AG][TC]|[CG][AT]{5}|[AT]{5}[CG]|[CG]{5}[AT]|[AT][CG]{5}'
  verboten_matches = re.findall(verb_strs + verb_pats, seq_str)
  return verboten_matches

def print_plots(file_name='James_Domains'):
  return print_th(file_name, plots=True)

def print_th_en(file_name='James_Domains'):
  return print_th(file_name, print_scores=True)

def get_th_scores(file_name='James_Domains'):
  return print_th(file_name)

def print_th(file_name='James_Domains', plots=False, print_scores=False):
  import matplotlib
  # matplotlib.use('QT4Agg') 
  import matplotlib.pyplot as plt
  
  import stickydesign as sd
  import energyfuncs_james as ej
  
  ef_j = ej.energyfuncs()
  ef_j.targetdG = 7.7
  
  toe_holds = sd.endarray(get_th_from_domains(file_name), 'TD')
  
  if plots:
    e_barplot(toe_holds, ef_j)
    plt.savefig('binding_energy.png')
    plt.close()

  e_vec = ef_j.matching_uniform(toe_holds)
  e_vec_ext = ef_j.th_external_dG(toe_holds)
  e_vec_int = ef_j.th_internal_dG(toe_holds)
  e_vec_avg = (e_vec_ext + e_vec_int) / 2
  e_vec_all = np.concatenate( (e_vec_int, e_vec_ext), 1)
  e_avg = e_vec_all.mean()
  e_rng = e_vec_all.max() - e_vec_all.min()
  
  if print_scores:
    output = 'TH Average Energy,TH Energy Range \n' + str(e_avg) + ',' + str(e_rng)
    f = open('TH_energies.txt', 'w')
    f.write(output)
    f.close()
  return [e_avg, e_rng]
