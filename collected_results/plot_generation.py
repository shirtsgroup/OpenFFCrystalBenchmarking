import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Script to generate plots of results

# import data
off1 = pd.read_csv('minimization_results_openFF1.csv')
off2 = pd.read_csv('minimization_results_openFF2.csv')
off1_box = pd.read_csv('minimization_results_box_minimized_OFF1.csv')
off2_box = pd.read_csv('minimization_results_box_minimized_OFF2.csv')

# convert RMSD values to numpy arrays
off1_rmsd = off1['RMSD'].str.strip('[]').to_numpy(dtype='float64')
off2_rmsd = off2['RMSD'].str.strip('[]').to_numpy(dtype='float64')
off1_box_rmsd = off1_box['RMSD'].str.strip('[]').to_numpy(dtype='float64')
off2_box_rmsd = off2_box['RMSD'].str.strip('[]').to_numpy(dtype='float64')

fig, axs = plt.subplots(2, 1)
axs[0].hist((off1_rmsd,off2_rmsd))
axs[0].set_xlabel('RMSD')
axs[0].set_ylabel('Frequency')
axs[0].legend(['OpenFF 1.0.0', 'OpenFF 2.0.0'])
axs[1].hist((off1_box_rmsd, off2_box_rmsd))
axs[1].set_xlabel('RMSD')
axs[1].set_ylabel('Frequency')
axs[1].legend(['OpenFF 1.0.0', 'OpenFF 2.0.0'])

fig.tight_layout(h_pad=3)
plt.savefig('RMSD_comparison.png')


plt.show()
