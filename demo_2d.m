%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  ACM TOG - J. Martinez Esturo, C. Rössl, and H. Theisel
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath util;

%% Compile mex file for polar decomposition
cd util; mex('-O', 'polardecomp.cpp'); cd ..;

%% Load mesh and show interactive 2D deformation window
beta = 0.33;

mesh = Mesh2D('meshes/giraffe2.off'); hidxs = [63,134,275,345,481,604,631,652];

def_asap   = Deform2DConformal(mesh,beta);
def_finite = Deform2DFinite   (mesh,beta,true);
def_arap   = Deform2DARAP     (mesh,beta,15);

viewctrl = DeformViewControl2D(def_arap);

viewctrl.addDeform(def_asap);
viewctrl.addDeform(def_finite);

viewctrl.show_energy = true;
viewctrl.show_energy_max = 6;

viewctrl.show();
viewctrl.initHandles(hidxs);
