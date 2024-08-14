# About this fork.

In this fork of openmm I am exploring the defintion of new types of forces. As of now I have created a CutoffAngleForce for the Reference platform. The idea behind it is that you can define a pair of atoms of which the distance is going to determine how strong the normal harmonic angle force is going to be. The cutoff function is sigmoidal to reduce artifacts during simulation. Notice that you can also use the cutoff value to entirely toggle the force, by setting it to 
a value <=0. 

I plan to implement a similar CutoffTorsionForce, and try to at least implement both for CUDA/OpenCL so it can be used in more intensive simulation.

The reason I am implementing these particular forces, is that although they can be easily implemented by following combination:

```python
distance = CustomBondForce('r')
angle_force = HarmonicAngleForce()
# add your angles
force = CustomCVForce('1/(1+exp(distance - cutoff))*angle_term')
force.addCollectiveVariable('distance', distance)
force.addCollectiveVariable('angle_term', angle_force)
force.addGlobalParameter('cutoff', 0.2*unit.nanometer)
```

But CustomCVForce is kind of a 'bubble' in memory and consumes lots of it (at least of shared memory on a gpu), so most of the time the system simulated like this has to be smaller than a standard MD simulation could be nowadays, or it has to be run on a CPU which is much slower. 


[![GH Actions Status](https://github.com/openmm/openmm/workflows/CI/badge.svg)](https://github.com/openmm/openmm/actions?query=branch%3Amaster+workflow%3ACI)
[![Conda](https://img.shields.io/conda/v/conda-forge/openmm.svg)](https://anaconda.org/conda-forge/openmm)
[![Anaconda Cloud Badge](https://anaconda.org/conda-forge/openmm/badges/downloads.svg)](https://anaconda.org/conda-forge/openmm)

## OpenMM: A High Performance Molecular Dynamics Library

Introduction
------------

[OpenMM](http://openmm.org) is a toolkit for molecular simulation. It can be used either as a stand-alone application for running simulations, or as a library you call from your own code. It
provides a combination of extreme flexibility (through custom forces and integrators), openness, and high performance (especially on recent GPUs) that make it truly unique among simulation codes.  

Getting Help
------------

Need Help? Check out the [documentation](http://docs.openmm.org/) and [discussion forums](https://simtk.org/forums/viewforum.php?f=161).
