%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Mesh < handle
  %MESH

  properties
    % user data
    data = []

    uv % parametrization assumed setable
  end

  properties (SetAccess = protected)
    nv, nt
    t

    uv0

    C  = []
    e  = []
    ie = [], iet = []
    bv = []

    A0, n0

    Mfull, Mbary, Mvoronoi
  end

  methods(Static)

  end

  methods
    function obj = Mesh(p,t,uv)
      assert(size(t,1) == 3);

      obj.t  = t;

      obj.nt = size(obj.t,2);

      obj.data = containers.Map;

      if isempty(uv),
        %use normalized projection
        uv = p(1:2,:);

        minu = min(uv(1,:));
        maxu = max(uv(1,:));
        uv(1,:) = (uv(1,:) - minu) ./ (maxu - minu);

        minv = min(uv(2,:));
        maxv = max(uv(2,:));
        uv(2,:) = (uv(2,:) - minv) ./ (maxv - minv);
      end

      %only use 2D uv's
      uv = uv(1:2,:);

      obj.uv  = uv;
      obj.uv0 = uv;
    end

    function C = get.C(obj)
      if(isempty(obj.C))
        obj.C = sparse(obj.t([1 2 3],:), obj.t([2 3 1],:), 1);
      end

      C = obj.C;
    end

    function e = get.e(obj)
      if(isempty(obj.e))
        obj.e = edges(obj.t, obj.C);
      end

      e = obj.e;
    end

    function ie = get.ie(obj)
      if(isempty(obj.ie))
        [obj.ie,obj.iet] = internaledges(obj.t, obj.C);
      end

      ie = obj.ie;
    end

    function iet = get.iet(obj)
      if(isempty(obj.iet))
        [obj.ie,obj.iet] = internaledges(obj.t, obj.C);
      end

      iet = obj.iet;
    end

    function bv = get.bv(obj)
      if(isempty(obj.bv))
        obj.bv = boundaryvertices(obj.t, obj.C);
      end

      bv = obj.bv;
    end

    function vidxs = closestVertices(obj,c)
      assert(size(c,1) == size(obj.p,1));

      nc = size(c,2);

      vidxs = zeros(1,nc);

      for i=1:nc,
        d = obj.p-repmat(c(:,i),1,obj.nv);
        d = sum(d.*d,1);

        [~,I] = min(d);
        vidxs(i) = I;
      end
    end

    function vidxs = closestParamVertices(obj,c)
      assert(size(c,1) == size(obj.uv,1));

      nc = size(c,2);

      vidxs = zeros(1,nc);

      for i=1:nc,
        d = obj.uv-repmat(c(:,i),1,obj.nv);
        d = sum(d.*d,1);

        [~,I] = min(d);
        vidxs(i) = I;
      end
    end

    function A0 = get.A0(obj)
      if(isempty(obj.A0)),
        [obj.n0, obj.A0] = trinormals(obj.p3_0, obj.t);
      end

      A0 = obj.A0;
    end

    function n0 = get.n0(obj)
      if(isempty(obj.n0)),
        [obj.n0, obj.A0] = trinormals(obj.p3_0, obj.t);
      end

      n0 = obj.n0;
    end

    function Mfull = get.Mfull(obj)
      if(isempty(obj.Mfull)),
        obj.Mfull = Mesh.massmatrix(obj.p3_0, obj.t, 'full', obj.A0);
      end

      Mfull = obj.Mfull;
    end

    function Mbary = get.Mbary(obj)
      if(isempty(obj.Mbary)),
        obj.Mbary = Mesh.massmatrix(obj.p3_0, obj.t, 'barycentric', obj.A0);
      end

      Mbary = obj.Mbary;
    end

    function h = drawUV(obj)
      h = trisurf(obj.t', ...
                  obj.uv(1,:), ...
                  obj.uv(2,:), ...
                  zeros(1,obj.nv), ...
                  zeros(1,obj.nv));
    end
  end

  methods (Abstract)
    h = draw(obj)
  end

  methods(Static)
    function M = massmatrix(p,t,type,A)
      %% This code is by Alex Jacobson (jacobson@inf.ethz.ch)

      %% Set up mesh mass matrix.
      %  - type: one of {'full', 'barycentric'}
      %  - A:    triangle areas (optional)

      if nargin<4
        [~, A] = trinormals(p, t);
      end

      dblA = 2.*A; %doubled area because we sum each triangle twice

      % renaming indices of vertices of triangles for convenience
      i1 = t(1,:); i2 = t(2,:); i3 = t(3,:);

      if strcmp(type,'full')
        % arrays for matrix assembly using 'sparse'
        % indices and values of the element mass matrix entries in the order
        % (1,2), (2,1),(2,3), (3,2), (3,1), (1,3) (1,1), (2,2), (3,3);
        i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
        j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
        offd_v = dblA/24.;
        diag_v = dblA/12.;
        v = [offd_v,offd_v, offd_v,offd_v, offd_v,offd_v, diag_v,diag_v,diag_v];
        M = sparse(i,j,v, size(p,2),size(p,2));
      elseif strcmp(type,'barycentric')
        % only diagonal elements
        i = [i1 i2 i3];
        j = [i1 i2 i3];
        diag_v = dblA/6.;
        v = [diag_v, diag_v, diag_v];
        M = sparse(i,j,v, size(p,2), size(p,2));
      else
        error('bad massmatrix type')
      end
    end
  end
end

