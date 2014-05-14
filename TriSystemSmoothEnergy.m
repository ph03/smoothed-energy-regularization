%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef TriSystemSmoothEnergy < handle
  %TRISYSTEMSMOOTHENERGY

  properties (SetAccess = protected)
    AAs
    bbs

    E = []
    b = []
  end

  properties (Access = protected)
    Wn   = []
    EtWn = []
    Um   = []
    vm   = []

    en   = []

    d    = []

    nv   = []
    nt   = []

    %for later total and local energy evaluation
    An   = []
    Bn   = []
    Dn   = []
  end

  methods
    function obj = TriSystemSmoothEnergy(mesh,beta, en, d,c)
      %% Generic Smoothed Energy Regularization
      %  This function implements the generic regularization of
      %
      %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %  Smoothed Quadratic Energies on Meshes
      %  J. Martinez Esturo, C. Rössl, and H. Theisel
      %
      %  ACM Transactions on Graphics 2014
      %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      %  For a given mesh and smoothing value the generic norm W_n is
      %  setup.
      %
      %  Using the functions updateLHS() and updateRHS() the actual
      %  discretized energy can be set up.
      %
      %  Its normal equation are given by the properties obj.AAs * u = obj.bbs.
      %
      %% Input:
      % mesh  - triangle mesh to apply energy smoothing to
      %         (mesh.p3_0 is used for integration domain / vertex coordinates)
      % beta  - smoothing factor in [0,1]
      % en    - number of quadratic error components per triangle
      % d     - dimension of unknown coefficients at vertices
      % c     - number of independent functions to setup system rhs
      %         (defaults to 1)
      %% Output:
      % AAs   - (property) the normal equations of the energy system
      % bbs   - (property) the rhs of the normal equations of the energy system

      if nargin < 6, c = 1; end

      %% number of error components
      obj.en = en;

      obj.d = d;

      obj.nv = mesh.nv;
      obj.nt = mesh.nt;

      %% compute triangle areas for replicated area weighting
      obj.An = sparse(1:(en*obj.nt), 1:(en*obj.nt), repmat(mesh.A0,en,1));

      ie  = mesh.ie;

      %% get internal edge lengths
      ie_lengths = mesh.p3_0(:,ie(2,:)) - mesh.p3_0(:,ie(1,:));
      ie_lengths = sqrt(sum(ie_lengths .* ie_lengths, 1));

      %% setup internal edge difference operator replicated over the en error
      %  components for each triangle
      nie = size(ie,2);

      obj.Dn = TriSystemSmoothEnergy.setupEdgeDifferenceOp(mesh,en);

      %% replicated edge length diagonal matrices
      obj.Bn = sparse(1:(en*nie), 1:(en*nie), repmat(ie_lengths,en,1));

      %% final smoothing inner product
      obj.Wn = (1-beta) .* obj.An + beta .* obj.Dn'*obj.Bn*obj.Dn;

      %% init rhs if b == 0
      obj.bbs = zeros(d *obj.nv,c);
      obj.b   = zeros(en*obj.nt,c);
    end

    function updateLHS(obj,E)
      assert(size(E,1) == obj.nt*obj.en);
      assert(size(E,2) == obj.nv*obj.d);

      obj.E = E;
      obj.EtWn = E'*obj.Wn;

      obj.AAs = obj.EtWn*E;
    end

    function updateRHS(obj,b)
      assert(size(obj.E,1) == obj.nt*obj.en);
      assert(size(obj.E,1) == size(b,1));

      obj.b   = b;
      obj.bbs = obj.EtWn*b;
    end

    function [terrEnergy,errEnergy,errSmoothn,errTotal] = ...
        evalTriangleErrors(obj,u)

      % local problem per triangle errors
      e = obj.E * u - obj.b;
      terrEnergy = sum(e.*e,2);
      terrEnergy = sum(reshape(terrEnergy,obj.en,[]));

      if nargout < 2, return; end;

      % integrated local errors
      errEnergy = trace(e'*obj.An*e);

      if nargout < 3, return; end;

      % integrated smoothness errors
      errSmoothn = obj.Dn*e;
      errSmoothn = trace(errSmoothn'*obj.Bn*errSmoothn);

      errTotal = errEnergy + errSmoothn;
    end
  end

  methods(Static, Access = private)
    function D = setupEdgeDifferenceOp(mesh,n)
      ie  = mesh.ie;
      iet = mesh.iet;

      nie = size(ie,2);

      D_i = [1:(n*nie),1:(n*nie)];
      D_j = repmat([n.*(iet(1,:)-1)+1,n.*(iet(2,:)-1)+1],n,1) + ...
              repmat((0:(n-1))',1,2*nie);
      D_v = [ones(n*size(iet,2),1),-1.*ones(n*size(iet,2),1)];
      D   = sparse(D_i, D_j, D_v, n * nie, n * mesh.nt);
    end
  end
end
