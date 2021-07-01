import numpy as np
from time import time
import sys
import os
import qdynos as qd

def quartic(x):
  return 0.01*x**4. - 0.5*x**2. + .1*x

def gauss(x,x0,sig,p0):
  return np.exp(-(x-x0)**2./2./sig**2.)*np.exp(-1.j*p0*x)

def basis(x,a,dx):
  return np.sqrt(dx)*np.sinc((x-a)/dx)

def norm(psi):
  return np.dot(psi.conj().T,psi)[0,0]

def make_model_stuff():

  hbar  = 1.0
  m   = 0.2
  kT  = 0.5

  xstar = 8.0
  dx  = 0.05
  x = np.arange(-xstar,xstar+dx,dx)
  interval = x[-1]-x[0]
  a = x[0]
  nstates = len(x)

  eta = 0.01
  ns = 10

  # Hamiltonian #
  H = np.zeros((nstates,nstates),dtype=complex)

  # potential
  for i in range(nstates):
    H[i,i] = quartic(x[i])

  # kinetic
  Tk = np.zeros((nstates,nstates))
  for i in range(nstates):
    for j in range(nstates):
      pre = hbar**2.*((-1.)**float(i-j))/2./m/dx**2.
      if i==j:
        Tk[i,j] = pre*np.pi**2./3.
      else:
        Tk[i,j] = pre*2./float(i-j)**2.
  H += Tk

  # position operator
  xop = np.diag(x)

  # diagonalize Hamiltonian and do transformations
  w,v = np.linalg.eigh(H)
  omega_c = w[2]-w[1]

  # add reorganization energy
  H += np.dot(xop,xop)*(eta*omega_c/np.pi)
  w,v = np.linalg.eigh(H)
  H = np.dot(v.conj().T,np.dot(H,v))
  xop = np.dot(v.conj().T,np.dot(xop,v))

  return hbar,kT,ns,H[:ns,:ns],xop[:ns,:ns],omega_c,eta

def main(argv):

  task = int(argv[1])
  tobs = int(argv[2])

  hbar,kT,ns,H,xop,omega_c,eta = make_model_stuff()

  # simulation stuff
  t_init = 0.
  dt = 1.0
  t_final = tobs*dt
  times0 = np.arange(t_init,t_final,dt)
  ntraj = int(argv[3])

  # expectation values
  e_ops = []
  e_ops.append( xop )

  # make system and bath
  bath = [qd.OhmicExp(eta, omega_c, kT, op=xop)]
  my_ham = qd.Hamiltonian(H, baths=bath)

  # make Lindblad class from Redfield theory
  gams , Ls , dynamics = qd.make_lindblad_from_rf(my_ham, sparse=False, adjoint=False)
  
  # run some setup things
  results = qd.Results(tobs=len(times0), store_states=True)
  options = qd.Options(unraveling=True)
  dynamics.dt = dt
  dynamics.setup(gams, Ls, options=options, results=results, eig=False, sparse=False, really_sparse=False)

  # window label
  lamda = int(argv[4])

  scratch_path = 'tobs%d/window%d/'%(tobs,lamda)
  try:
    os.mkdir('tobs%d'%(tobs))
  except:
    print('tobs dir made already')
  try:
    os.mkdir('tobs%d/window%d'%(tobs,lamda))
  except:
    print('window dir made already')

  vib = np.zeros((ns,ns))
  for i in range(ns):
    vib[i,i] = float(i)

  # projector defining h_A indicator function
  projector_a = np.zeros((ns,ns))
  projector_a[0,0] = 1.
  projector_a[2,2] = 1.

  # projector defining h_B indicator function based on the window (lamda)
  projector_b = None
  if lamda == 0:
    projector_b = np.zeros((ns,ns))
    projector_b[2,2] = 1.
    projector_b[4,4] = 1.
  elif lamda == 1:
    projector_b = np.zeros((ns,ns))
    for i in range(4,ns):
      projector_b[i,i] = 1.
  elif lamda == 2:
    projector_b = np.zeros((ns,ns))
    projector_b[5,5] = 1.
    projector_b[3,3] = 1.
  elif lamda == 3:
    projector_b = np.zeros((ns,ns))
    projector_b[3,3] = 1.
    projector_b[1,1] = 1.
  elif lamda == 4:
    projector_b = np.zeros((ns,ns))
    projector_b[1,1] = 1.
  elif lamda == 5:
    projector_b = np.zeros((ns,ns))
    projector_b[0,0] = 1.
    projector_b[2,2] = 1.

  # read in initial trajectory #
  psi0 = []
  fpsi = open(argv[5],'r')
  for i in range(tobs):
    psi = np.zeros((ns,1),dtype=complex)
    line = fpsi.readline().split()
    for j in range(ns):
      psi[j,0] = complex(line[j])
    psi0.append( psi.copy() )
  fpsi.close()

  # observables of current trajectory
  ha_prev = np.around(((np.dot(psi0[0].conj().T,np.dot(projector_a,psi0[0]))[0,0]).real),6)
  hb_prev = np.around(((np.dot(psi0[-1].conj().T,np.dot(projector_b,psi0[-1]))[0,0]).real),6)
  viba_prev = np.around(((np.dot(psi0[0].conj().T,np.dot(vib,psi0[0]))[0,0]).real),6)
  vibb_prev = np.around(((np.dot(psi0[-1].conj().T,np.dot(vib,psi0[-1]))[0,0]).real),6)

  # random number generator to choose point in trajectory to shoot
  random_shooter = np.random.RandomState(int(time()+100*task))

  # random number generator for dynamics
  rng = np.random.RandomState(int(time()+100*task+np.random.randint(100000)))

  ftshoot = open(scratch_path+'t_shoots_tobs_%d_task_%d.dat'%(tobs,task),'w')
  fhb = open(scratch_path+'hb_list_tobs_%d_task_%d.dat'%(tobs,task),'w')

  # status updater
  every = int(ntraj/10)

  # statistics counters of the number of accepted trajectories forward and backward
  nbaccept = 0
  nfaccept = 0

  # this is akin to the number of monte carlo sweeps, we'll generate technically more trajectories, but we only
  # look at the resulting trajectory after all initial positions have been sampled over
  for i in range(ntraj):

    print(i,nbaccept,nfaccept)
    if i%every==0:
      print("%.0f percent done"%(int(10*(i/every))))
      fpsi = open(scratch_path+argv[6],'w')
      for j in range(tobs):
        for k in range(ns):
          fpsi.write(str(psi0[j][k,0])+' ')
        fpsi.write('\n')
      fpsi.close()

    # sample over all the times along the trajectory
    for j in range(tobs):
      # pick random point along trajectory to backward shoot from #
      t_shoot = random_shooter.randint(2,tobs)
      times = times0[:t_shoot]
      output = dynamics.integrate_trajectory(psi0[t_shoot].copy(), times, rng)

      # compute endpoint observables #
      ha_new = int(np.around(((np.dot(output.states[-1].conj().T,np.dot(projector_a,output.states[-1]))[0,0]).real)))
      vib_new = int(np.around(((np.dot(output.states[-1].conj().T,np.dot(vib,output.states[-1]))[0,0]).real)))

      # choose to accept or reject the backward shot
      if ha_new > 0.5:
        viba_new = vib_new
        ftshoot.write('%d b\n'%(t_shoot))
        ftshoot.flush()
        # accept
        nbaccept += 1
        ha_prev = ha_new
        viba_prev = vib_new
        # update trajectory
        for k in range(t_shoot):
          psi0[k] = output.states[-k].copy()

      # pick random point along trajectory to forward shoot from #
      t_shoot = random_shooter.randint(0,tobs-1)
      times = times0[t_shoot:]
      output = dynamics.integrate_trajectory(psi0[t_shoot].copy(), times, rng)

      # compute enpoint observables #
      hb_new = int(np.around(((np.dot(output.states[-1].conj().T,np.dot(projector_b,output.states[-1]))[0,0]).real)))
      vib_new = int(np.around(((np.dot(output.states[-1].conj().T,np.dot(vib,output.states[-1]))[0,0]).real)))

      if hb_new > 0.5:
        #if int(round(hb_new)) == 0:
        #  print(hb_new,vib_new)
        ftshoot.write('%d\n'%(t_shoot))
        ftshoot.flush()
        # accept trajectory
        nfaccept += 1
        hb_prev = hb_new
        vibb_prev = vib_new
        for k in range(t_shoot,tobs):
          psi0[k] = output.states[k-t_shoot].copy()

    fhb.write('%d %d %d %d\n'%(int(round(ha_prev)),int(round(viba_prev)),int(round(hb_prev)),int(round(vibb_prev))))
    fhb.flush()

  tot_traj = ntraj*tobs
  print(int(tot_traj),nbaccept,nfaccept)
  ftshoot.close()
  fhb.close()
  
  # print final traj
  fpsi = open(scratch_path+argv[6],'w')
  for j in range(tobs):
    for k in range(ns):
      fpsi.write(str(psi0[j][k,0])+' ')
    fpsi.write('\n')
  fpsi.close()

if __name__ == "__main__":
  main(sys.argv)
