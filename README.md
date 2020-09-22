# QCAT

A Fortran program of **Q**uasi-**C**lassical **A**diabatic **T**rajectory for collision of two molecules. (under development)

- **Author** : Wenbin, FAN (fanwenbin@shu.edu.cn)
- **Supervisor** : Yongle, LI

# Intro

* **Unit** : atomic unit (in code)

* **Language** : Fortran 95
* **Requirement** : Intel Fortran Compiler, MKL

# Usage

1. Clone all codes. 
2. Connect your PES in `Interface.f90`. 
3. Modify input file with any name. 
4. Modify compiler and input file name in the `Makefile`. 
5. `make`. 

Note that the test case is OHCH4 system, whose PES could be obtained from Prof. Bin JIANG upon request. [ J. Li and H. Guo, J Chem Phys **143**, 221103 (2015). ]

# Keyword in input file

* `temperature` real value (K), default 298.15 K. 
* `ensemble` *NVE*, *NVT*, Nose-Hoover thermostat
* `timestep` real value (fs), default 0.01 fs. 
* `ntraj` the number of trajectories, default 1. 
* `nstep` the number of steps in each trajectories. 
* `nprint` print frequency of system properties. 
* `nxyz` print frequency of trajectory. 

# Developing Progress

**Finished : **

1. **initialize molecule**
   1. remove center of mass
   2. rotate molecule to principal axis
2. **initialize trajectory**
   1. random rotation
   2. random thermal momentum
   3. random rotational state
   4. random vibrational state
   5. remove translational and angular momentum
3. **propagate trajectory**
   1. Velocity Verlet
   2. Nose-Hoover thermostat
4. **interface and utilities**
   1. interface to standard PES (coordinates and energy)
   2. gradient and hessian matrix, vibrational analysis
   3. pseudo-inversion of matrix

**TODO : **

1. the correctness of all kinds of momentum
2. analysis of final states and trajectory
3. treatment of ZPE