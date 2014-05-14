Smoothed Quadratic Energies on Meshes
=====================================
ACM Transactions on Graphics 2014 - J. Martinez Esturo, C. Rössl, and H. Theisel
--------------------------------------------------------------------------------

This repository contains a reference implementation of the *smoothed energy
regularizer for quadratic energies on triangular meshes* proposed in [Paper].

MATLAB Demo Implementation
--------------------------
The energy regularization reference implementation is given in
`TriSystemSmoothEnergy.m`.

As an example application, we provide a demo of interactive 2D deformations the
minimize the ARAP, lARAP, and ASAP energies.

Usage
-----
Clone the reposityory with
`git clone --recursive https://github.com/ph03/smoothed-energy-regularization`.

Run the `demo_2d.m` file in matlab to start.

Requirements
------------
A recent Matlab version with object-oriented language features is required. We
have tested the code on both Windows and Linux machines with Matlab versions
r2012a and r2013a.

Credits
-------
We link the sources of the Eigen library (see [Eigen]) as a submodule to
efficiently compute polardecompositions of deformation gradients.

The code makes use of Matlab packages by Alec Jacobson (see [IGLlab]) for mesh
loading. Please contact Alec at jacobson@inf.ethz.ch and Janick at
janick@mpi-inf.mpg.de, before using this code outside of an informal
setting, i.e., for comparisons in an academic paper.

All files copyright Janick Martinez Esturo 2014 (MIT License) unless
otherwise noted.

[Paper]:  mpi-inf.mpg.de/~jmartine/pubdetails/MartinezEsturo2014.html
[Eigen]:  eigen.tuxfamily.org
[IGLlab]: igl.ethz.ch
