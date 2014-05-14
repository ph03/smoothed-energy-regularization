%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Deform2DARAP < Deform2D
  %DEFORM2DARAP

  properties
    se = []
    constrSolver = []

    def_finite = []

    itmax = []
  end

  methods
    function obj = Deform2DARAP(mesh, beta, itmax)
      %% Constructor
      obj = obj@Deform2D(mesh, beta);

      obj.itmax = itmax;

      obj.name = 'ARAP';
    end

    function init(obj, hidxs)
      %% Setup finite deformation operator (some redundant computations atm...)
      mesh = obj.mesh;

      obj.def_finite = Deform2DFinite(mesh, obj.beta, false);
      obj.def_finite.csolverf = obj.csolverf;
      obj.def_finite.init(hidxs);

      %% Initialize local / global iteration

      % global smoothed energy is just a poisson system
      en = size(mesh.GP,1) / mesh.nt;

      obj.se = obj.smootherf(mesh,obj.beta, en, 1,1);
      obj.se.updateLHS(mesh.GP);

      %% Setup global solver

      % constrained handle indices
      obj.hidxs = hidxs;

      obj.constrSolver = obj.csolverf(obj.se.AAs,obj.hidxs,[]);
    end

    function [converged,u] = deform(obj, hcoords, itn)
      converged = true; u = [];
      if(isempty(obj.hidxs)), return; end

      if itn == 1,
        %% initialize with linearized finite deformation
        obj.def_finite.deform(hcoords, itn);
        converged = false;
        u=[];
      else
        %% local / global iteration
        mesh = obj.mesh;

        %% local: find best rotation per triangle deformation gradient
        x = reshape(mesh.p,[],1);

        F = mesh.GGP*x;

        if false,
          % matlab svd
          R = zeros(2,2*obj.mesh.nt);

          for t=1:mesh.nt,
            [U,S,V] = svd(reshape(F((4*(t-1)+1):(4*t),1),2,2));

            L = U*V';

            if S(1,1)*S(2,2) < 0,
              L = -1.*L;
            end

            R(:,(2*(t-1)+1):(2*t)) = L;
          end
        else
          % mex svd
          F = reshape(F,2,[]);

          opt.fast2D = false;
          R = polardecomp(F, opt);
        end

        %% global: fit deformation to these gradients

        % reshape rotations to columnwise blocks
        R = blockreshape(R',2,1,mesh.nt)';

        obj.se.updateRHS(R);
        u = obj.constrSolver.solve(hcoords', obj.se.bbs);

        mesh.p = u';

        if itn == obj.itmax,
          converged = true;
        else
          converged = false;
        end
      end
    end
  end
end
