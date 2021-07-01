# qtps
Example of performing a quantum transition path sampling calculation using [qdynos](https://github.com/addschile/qdynos).

There are two folders where you'll find relevant parts of an example calculation taken from J. Chem. Phys. 149, 214109 (2018).

- qtps - folder containing the code for running the calculation
- data - this folder contains example data for Fig. 9 from the above reference

# Running a qtps calculation

You'll first need to install [qdynos](https://github.com/addschile/qdynos), which is a little python package that I've written to do various quantum dynamics calculations. Feel free to smash that like button if you find the code interesting or useful!

As part of the code you'll find three files:

- quartic_model.py - This is the main driver code for the specific model. In this, you'll find code that creates the Hamiltonian using a sinc function DVR, code that sets up the Lindblad dynamics from a microscopic model (via the frequency projection of the dressed Redfield operators using the qdynos code), and the main TPS loop by which new trajectories are generated from old ones with Monte Carlo. I've tried to comment up the code significantly so that at each step it is clear what is going down. Something that I should mention is that because everything is done in the eigenbasis of the Hamiltonian and since the Lindblad operators are time-reversal invariant, you can do both forward and backward trajectories. If you have a more complicated initial condition or a more complicated completely positive trajectory equation (e.g. HOPS) you would have to be more careful about doing both forward and backward trajectories and might only be able to do forward trajectories, which is not a problem in this example since trajectories are very much instanton-like. 

- run.sh - you can run this with 'bash run.sh'. It contains an example of the command line that will run the tps calculation for a specific window. For the calculation in the paper, I ran this script changing the "task" variable in the command line from 0-9. I used ntraj = 15000, which generates 15000 trajectories from monte carlo sweeps and I used the final 10000 trajectories from each block for each window. The data generated and analysis of that data is described below.

- psi_init.txt - example of an initial trajectory for the tps calculation. It is always wise to use different initial trajectories for different block averaging, but as long as you're using different random number seeds, a long "equilibration" period (in the number of monte carlo sweeps) should take care of that. Also, when using different windows it's wise to make different initial trajectories that end within the window you want, though again with a long enough equilibration period that shouldn't be an issue.

Feel free to reach out if you have any questions, I'm more than happy to help!
