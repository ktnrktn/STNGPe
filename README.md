# STNGPe
simulation of the STN-GPe network in the basal ganglia

### Outline
This program conducted numerical simulation of sub thalamic nucleus (STN) and external globus pallidus (GPe) network in the basal ganglia.
The network consisted detailed neuron models of STN neurons, prototypic GPe neurons, and arkypallidal neurons.
The code is parallelized using MPICH.

---

### Requirement
* MPICH-3.3.2 (or later)
* C-compiler(gcc etc.)
* Perl (optional)

---

### Usage
The simulation consisted two phases.
Firstly, it is necessary to make configuration files for network connections between cell populations.
Then run network simulation program that reads the network configuration files.

While “netconn” is used in the first step, “stngpMPI” in the second step.

Compile with,
type
`% make netconn`
and 
`% make stngp`

First, run “netconn”
`% ./netconn arg1 … arg9`
To see what the arguments mean,
type without any arguments
`%./netconn`
To run a network simulation,
type
`% mpirun -np 8 ./stngpMPI arg1 … arg33`
Because the parameters are too much, use the Perl script “exec.pl” in which all the parameters are set.
`% ./exec.pl`

"neuronXXX.conf" files are configuration files for neuron models.
The file contains maximal conducance of ion channels.
"stngpMPI" reads these files and set the conductance parameters.
