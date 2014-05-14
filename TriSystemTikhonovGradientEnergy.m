%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef TriSystemTikhonovGradientEnergy < handle
  %TRISYSTEMTIKHONOVGRADIENTENERGY

  properties (SetAccess = protected)
    AAs
    bbs

    E = []
    b = []
  end

  properties (Access = protected)
    en   = []

    d    = []

    nv   = []
    nt   = []

    beta = []

    %for later total and local energy evaluation
    An   = []

    AT   = []
    T    = []
  end

  methods
    function obj = TriSystemTikhonovGradientEnergy(mesh,beta, en, d,c)
      %%
      % mesh  - triangle mesh to apply energy smoothing to
      %         (mesh.p3_0 is used for integration domain / vertex coordinates)
      % beta  - smoothing factor in [0,1]
      % en    - number of quadratic error components per triangle
      % d     - dimension of unknown coefficients at vertices
      % c     - number of independent functions to setup system rhs
      %         (defaults to 1)

      if nargin < 6, c = 1; end

      obj.beta = beta;

      %% number of error components
      obj.en = en;

      obj.d = d;

      obj.nv = mesh.nv;
      obj.nt = mesh.nt;

      %% compute triangle areas for replicated area weighting
      obj.An = sparse(1:(en*obj.nt), 1:(en*obj.nt), repmat(mesh.A0,en,1));

      %% regularization by gradient operator
      %  (we use 3D gradient vectors)
      A3 = sparse(1:(3*obj.nt), 1:(3*obj.nt), repmat(mesh.A0,3,1));
      G  = mesh.GP_3;

      perm = reshape(1:(d*obj.nv),d,[])';
      P    = sparse (1:(d*obj.nv),perm,1);

      for i=1:d,
        Tc{i}  = G;
        ATc{i} = A3;
      end

      obj.T  = spblkdiag(Tc) * P;
      obj.AT = spblkdiag(ATc);

      %% init rhs if b == 0
      obj.bbs = zeros(d*obj.nv,c);
      obj.b   = zeros(en*obj.nt,c);
    end

    function updateLHS(obj,E)
      assert(size(E,1) == obj.nt*obj.en);
      assert(size(E,2) == obj.nv*obj.d);
      obj.E = E;

      obj.AAs = (1-obj.beta) * obj.E'*obj.An*obj.E + ...
                   obj.beta  * obj.T'*obj.AT*obj.T;
    end

    function updateRHS(obj,b)
      assert(size(obj.E,1) == obj.nt*obj.en);
      assert(size(obj.E,1) == size(b,1));

      obj.b   = b;
      obj.bbs = (1-obj.beta) * obj.E'*obj.An*b;
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
      errSmoothn = obj.T*u;
      errSmoothn = trace(errSmoothn'*obj.AT*errSmoothn);

      errTotal = errEnergy + errSmoothn;
    end
  end
end
