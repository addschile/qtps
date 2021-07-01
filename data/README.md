This folder contains example data for Fig. 9 from J. Chem. Phys. 149, 214109 (2018).

The TPS calculation is repeated for a number of different trajectory lengths (tobs, the number of timesteps) and runs for a particular trajectory length is in folder tobs. 

Within each tobs folder there are a number of folders named 'windowX' where X is 0 to 5. Each window represents the window covered by the h_B operator and they should overlap adjacent windows (e.g. window 2 overlaps with windows 1 and 3).

Within each window are a bunch of files 'hb_list_tobs_X_task_Y.dat'. X is the tobs (200, 300, etc) and Y is an integer between 0 and 19 (inclusive) labelling an independent run of the TPS calculation (this was done to take advantage of supercomputing resources). Each file has lines of data of the form:

h_A(0) vib_state(0) h_B(tobs) vib_state(tobs)

h_A(0) is the expectation value of the h_A operator at time 0 of the trajectory (should always be 1)
h_B(tobs) s the expectation value of the h_B operator at time tobs of the trajectory (should always be 1)
vib_state(0) and vib_state(tobs) are the vibrational states at time 0 and tobs, respectively (vib_state(0) should always be 0 or 2, which are the vibrational states associated with the A state and vib_state(tobs) should be whatever states are in that particular window, e.g. for window 0 they are 2 and 4)

The data from all the different tasks (I only used 10 tasks, 0-9) includes "equilibration" and production, so the program 'make_data_files.py' is used to take the final portion of data and put into a single file called 'windowX.dat' that has the same form as mentioned above.

This vib_state(tobs) data is the data that is ultimately reweighted since that is what is biased in the TPS calculation. The reweighting is done in the script 'reweight_dists.py'. This code reads in all the data for a particular tobs and computes the probability densities of ending in a particular vibrational state, which is a biased distribution. These distributions are reweighted by hand by minimizing the difference between the distribution at the overlapping states in each window and then, once reweighted the new distributions are printed to a file 'reweighted_hist_X.dat' where X is a particular block of data that was used in this calculation to make block averages for the rate. The reweighted distribution data is printed out of the form:

vib_state(tobs) -ln(probability)

You'll notice some of the vib_state(tobs) are repeated, which is because the windows are overlapping in the vibrational state and the log of the probabilities should be equal for repeated vib_states(tobs)

The rates are finally calculated in 'compute_rates.py' as the ratio of the path partition functions, which is a log ratio of the reweighted probabilities. This ratio is computed for each window and printed out to 'path_partition_function_X.dat' (X is tobs) and then the average and standard deviation are printed at the bottom. Note the unit system the calculations are done with are different what was is reported in the paper. 

These final data points are organized for each tobs in the file 'path_ensemble_ratios.dat', which is plotted in 'plot_ensemble_ratio.py'.

For tobs = 200, the ratio of partition functions is printed out for each window in 'avg_reweighted_hist.dat' as a function of the average position of the vibrational wavefunction and plotted with 'make_fig2.py'. This is the final figure in the paper (Fig. 9).
