Smoothed Quadratic Energies on Meshes
=====================================

ACM TOG - J. Martinez Esturo, C. Rössl, and H. Theisel
------------------------------------------------------

MATLAB Demo Implementation by J. Martinez Esturo
------------------------------------------------
This archive contains a reference implementation of the smoothed energy
regularizer for triangular meshes (TriSystemSmoothEnergy.m).

As an example application we provide a demo of interactive 2D deformations the
minimize the ARAP, lARAP, and ASAP energies.

Please run the *demo_2d.m* file to start.

Requirements
------------
A recent Matlab version with object-oriented language features is required. We
have tested the code on both Windows and Linux machines with Matlab versions
r2012a and r2013a.

Credits
-------
We provide the sources of the Eigen library (see [Eigen]) used to
compute polardecompositions of deformation gradients efficiently.

The code makes use of Matlab packages by Alec Jacobson (see [IGLlab]) for mesh
loading. Please contact Alec at jacobson@inf.ethz.ch and Janick at
janick@mpi-inf.mpg.de, before using this code outside of an informal
setting, i.e., for comparisons in an academic paper.

All files copyright Janick Martinez Esturo 2014 unless otherwise noted.

[Eigen]:  eigen.tuxfamily.org
[IGLlab]: igl.ethz.ch

