## CHEME 7770 Example: Three gene memory motif
This repository contains the model code for a three-gene memory motif example implemented in the [Julia](http://julialang.org) programming language.

### Installation and Requirements
You can download this repository as a zip file, or `clone`/`pull` it by using the command (from the command-line):

	$ git pull https://github.com/varnerlab/CHEME7770-ThreeGeneModel-Example.git

or

	$ git clone https://github.com/varnerlab/CHEME7770-ThreeGeneModel-Example.git

The model code was machine generated using the [Gene Regulatory Network in Julia (JuGRN)](https://github.com/varnerlab/JuGRN-Generator) code generation system. The model code uses several [Julia](http://julialang.org) packages:

Package | Description | Command
--- | --- | ---
ODE | Contains the ``ode23`` subroutine to solve the model equations | Pkg.add("ODE")
PyPlot | Used to make figures (assume you have Python installed) | Pkg.add("PyPlot")

### Solve the model equations?
The model equations can be solved by executing one of the predefined driver routines.
``Washout.jl`` solves the model equations where we add, and then remove the inducer for the default network,
while ``WashoutKO.jl`` solves the same case with a network in which the positive feedback between ``P2`` and ``P3`` is broken.
The model equation, and simulation code is are contained in the [src](https://github.com/varnerlab/CHEME7770-ThreeGeneModel-Example/tree/master/src) directory. 
