# TOF

This Python-2 script is used to run 1D time-of-flight simulations to compare different types of energy discretization, including the: 

* multigroup (MG) method
* finite-element-with-discontiguous-support (FEDS) method
* finite-element-with-discontiguous-support method where each subelement is given a separate velocity (FEDS-sub)

This script was written to solve a single 1D analytical test problem.

## 1D analytical time-of-flight problem ##

This neutron tranport problem includes:
* a slab of U-235 for <img src="http://latex.codecogs.com/gif.latex?x\in[0,A)" border="0"/>
* a slab of natural iron for <img src="http://latex.codecogs.com/gif.latex?x\in[A,A+B)" border="0"/>
* a vacuum for <img src="http://latex.codecogs.com/gif.latex?x\in[A+B,D)" border="0"/>
* a point detector at <img src="http://latex.codecogs.com/gif.latex?x=D" border="0"/>

where A = 0.05 cm, B = 2.06 cm, and D = 100 cm. The slab of U-235 has a fission source uniformly-distributed in x, where all neutrons are emitted in the forward <img src="http://latex.codecogs.com/gif.latex?\mu=1" border="0"/> direction. Both the U-235 and the iron slabs are pure absorbers, and the total macroscopic cross section is used to calculate the attenuation for both slabs. 

We assume the fission source <img src="http://latex.codecogs.com/gif.latex?\chi_{235}(E)" border="0"/> is only pulsed for times <img src="http://latex.codecogs.com/gif.latex?t\in[0,\tau]" border="0"/>, and measure the time-dependent detector response.

--------------

More details of this code are explained in the documentation folder.
