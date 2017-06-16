import sys
import os

Un = 4e-04
Ub = 7.5e-05
s0 = 7e-03
Xc = 5e-02
g = 1.0/2/Xc

num_replicates = 108

sys.stderr.write("Running forward-time simulations!\n")

os.system("additional_data/macroscopic_epistasis_simulations/code/simulate_wiser_timecourse %d %g %g %g %g > additional_data/macroscopic_epistasis_simulations/wiser_timecourse_output.txt" % (num_replicates, Ub, s0, g, Un))

sys.stderr.write("Done!\n")