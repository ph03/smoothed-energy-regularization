%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Deform2DFinite < Deform2D
  %DEFOR2DMFINITE

  properties
    se = []
    constrSolver = []

    mode = [];
  end

  methods
    function obj = Deform2DFinite(mesh, beta, mode_lrigid_or_conf)
      %% Constructor
      obj = obj@Deform2D(mesh, beta);
      obj.mode = mode_lrigid_or_conf;

      if mode_lrigid_or_conf,
        obj.name = 'lARAP';
      else
        obj.name = 'ASAP';
      end
    end

    function init(obj, hidxs)
      %% Setup element permutation opertor P
      mesh = obj.mesh;

      %% Setup material behaviour term

      if ~obj.mode,
        %% conformal <=> as-similar-as-possible
        Mwq = eye(4) - 0.5 .* [1, 0, 0,1;
                               0, 1,-1,0;
                               0,-1, 1,0;
                               1, 0, 0,1];

        Mq = blockfill(4,4,mesh.nt,repmat(Mwq,1,mesh.nt));

        %% Setup integrated / smoothed system
        E  = Mq*mesh.GGP;

        en = size(E,1) / mesh.nt;

        obj.se = obj.smootherf(mesh,obj.beta, en, 2,1);
        obj.se.updateLHS(E);
      else
        %% linearized rotations
        Mwq = 1.0 .* [1,  0,  0,0;
                      0,1/2,1/2,0;
                      0,1/2,1/2,0;
                      0,  0,  0,1];
        Mwr = [1;0;0;1];

        Mq = blockfill(4,4,mesh.nt,repmat(Mwq,1,mesh.nt));

        %% Setup integrated / smoothed system
        E  = Mq*mesh.GGP;
        b  = repmat(Mwr,mesh.nt,1);

        en = size(E,1) / mesh.nt;

        obj.se = obj.smootherf(mesh,obj.beta, en, 2,1);
        obj.se.updateLHS(E);
        obj.se.updateRHS(b);
      end

      %% Setup solver

      % constrained handle indices
      obj.hidxs = hidxs;
      cidxs = reshape([(hidxs-1)*2+1; (hidxs-1)*2+2],1,[]);

      obj.constrSolver = obj.csolverf(obj.se.AAs,cidxs,obj.se.bbs);
    end

    function [converged,u] = deform(obj, hcoords, ~)
      converged = true; u = [];
      if(isempty(obj.hidxs)), return; end

      u = obj.constrSolver.solve(reshape(hcoords,[],1));

      obj.mesh.p = reshape(u,2,[]);
    end
  end
end
