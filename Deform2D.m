%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  ACM TOG - J. Martinez Esturo, C. Rössl, and H. Theisel
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Deform2D < handle
  %DEFORM2D

  properties
    % smoothing weight
    beta

    % solver creation lambda function
    csolverf

    % smoother creation lambda function
    smootherf
  end

  properties (SetAccess = protected)
    % mesh to deform
    mesh

    % handle indices
    hidxs = []

    % name of the deformation method
    name
  end

  methods
    function obj = Deform2D(mesh, beta)
      %% Constructor
      obj.mesh = mesh;
      obj.beta = beta;

      %% default solver
      obj.csolverf = @(A,cidxs,b)CholConstrSolver(A,cidxs,b);

      %% default smootherf
      obj.smootherf = @(mesh,beta, en, d,c) ...
        TriSystemSmoothEnergy(mesh,beta, en, d,c);
    end

    function set.beta(obj, beta)
      obj.beta = beta;

      %% re-init deform implementation
      if(~isempty(obj.hidxs))
        obj.init(obj.hidxs);
      end
    end

    function set.csolverf(obj, csolverf)
      obj.csolverf = csolverf;

      %% re-init deform implementation
      if(~isempty(obj.hidxs))
        obj.init(obj.hidxs);
      end
    end

    function set.smootherf(obj, smootherf)
      obj.smootherf = smootherf;

      %% re-init deform implementation
      if(~isempty(obj.hidxs))
        obj.init(obj.hidxs);
      end
    end
  end

  methods (Abstract)
    init  (obj, hidxs)
    [converged,u] = deform(obj, hcoords, itn)
  end
end
