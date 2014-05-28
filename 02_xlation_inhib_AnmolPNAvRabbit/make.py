import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

tbl = pd.read_excel('in/fluor_rab_crassa.xlsx','results')
avg_rab = np.mean(tbl[['rab1','rab2', 'rab3']],1)
std_rab = np.std(tbl[['rab1','rab2', 'rab3']],1)
avg_crs = np.mean(tbl[['cras1','cras2']],1)
std_crs = np.std(tbl[['cras1','cras2']],1)

exp_mean = (avg_rab[0], avg_crs[0])
ctrl_mean = (avg_rab[1], avg_crs[1])
exp_std = (std_rab[0], std_crs[0])
ctrl_std = (std_rab[1], std_crs[1])

fig = plt.figure()
ax = fig.add_subplot(111)
x = np.arange(2)
width = 0.35
bar_exp = ax.bar(x, exp_mean, width, color='r', yerr=exp_std)
bar_ctrl = ax.bar(x+width, ctrl_mean, width, color='y', yerr=ctrl_std)
ax.set_ylabel('Fluorescence (a.u.)')
ax.set_title('Ribosomal Inhibition in Rabbit, Crassa')
ax.set_xticks(x+width)
ax.set_xticklabels(('Rabbit','Crassa'))
ax.legend((bar_exp,bar_ctrl), ('2uM PNA','0uM PNA'))

plt.savefig('out/02_pna_v_rabcras.pdf')