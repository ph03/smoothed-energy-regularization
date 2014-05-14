%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Mesh2D < Mesh
  %MESH2D

  properties
    p % positions assumed setable
  end

  properties (SetAccess = protected)
    p_0

    p2, p3
    p2_0, p3_0

    % (deformation) gradient operators
    GP
    GGP
    GP_3

    texture = []
    force_texture = false
  end

  methods(Static)
  end

  methods
    function obj = Mesh2D(p,t,uv)
      %% Constructor
      if nargin < 2,
        [p,t,uv] = readMesh(p);
      end

      obj = obj@Mesh(p,t,uv);

      assert(size(p,1) == 3 || size(p,1) == 2);

      obj.p   = p(1:2,:);
      obj.p_0 = obj.p;

      obj.nv = size(obj.p,2);
    end

    function p2   = get.p2(obj),     p2 =  obj.p;   end
    function p3   = get.p3(obj),     p3 = [obj.p;   zeros(1,obj.nv)]; end
    function p2_0 = get.p2_0(obj), p2_0 =  obj.p_0; end
    function p3_0 = get.p3_0(obj), p3_0 = [obj.p_0; zeros(1,obj.nv)]; end

    function GGP = get.GGP(obj)
      if(isempty(obj.GGP))
        obj.updateGradOps2D();
      end

      GGP = obj.GGP;
    end

    function GP = get.GP(obj)
      if(isempty(obj.GP))
        obj.updateGradOps2D();
      end

      GP = obj.GP;
    end

    function GP_3 = get.GP_3(obj)
      if(isempty(obj.GP_3))
        GP_3 = Mesh3D.getGradOps3D(obj);

        obj.GP_3  = GP_3;
      end

      GP_3 = obj.GP_3;
    end

    function setTexture(obj,texturefname)
      if isempty(texturefname),
        obj.texture = [];
        return;
      end

      obj.texture = imread(texturefname);

      % setup uv coordinates
      uv = [obj.p(1,:);obj.p(2,:)];

      minx = min(uv(1,:)); maxx = max(uv(1,:));
      miny = min(uv(2,:)); maxy = max(uv(2,:));

      uv(1,:) = (uv(1,:) - minx) ./ (maxx - minx);
      uv(2,:) = (uv(2,:) - miny) ./ (maxy - miny);

      obj.uv = [1-uv(2,:);uv(1,:)];
    end

    function setForceTexture(obj,force)
      obj.force_texture = force;
    end

    function h = draw(obj)
      textured = obj.force_texture;

      if isempty(obj.texture)
        textured = false;
      end

      if ~textured,
        h = trisurf(obj.t', ...
                    obj.p(1,:), ...
                    obj.p(2,:), ...
                    zeros(1,obj.nv));
      else
        % textures suck hard in MATLAB :(
        h = patcht(obj.t',[obj.p;zeros(1,obj.nv)]',obj.t',obj.uv', ...
                   obj.texture);
      end
    end
  end

  methods (Access = protected)
    function updateGradOps2D(obj)
      %% Compute combined deformation gradient operator GGP

      % Setup per element deformation gradient operator
      Gvalues  = zeros(2,3*obj.nt);
      GGvalues = zeros(4,6*obj.nt);

      for t=1:obj.nt,
        % solve 2x2 system for triangle gradient operator
        %   and permute into solution

        p_0 = obj.p_0(:,obj.t(:,t));

        % 2d version of Botsch et al. - Deformation
        %   Transfer for Detail-Preserving Surface Editing
        % -> triangle gradient operator
        G = [(p_0(:,2) - p_0(:,1))'; (p_0(:,3) - p_0(:,1))'] \ [[-1;-1],eye(2)];

        Gvalues(:,(3*(t-1)+1):(3*t)) = G;

        % deformation gradient operator for blocked coordinates
        GG = [G,zeros(2,3);zeros(2,3),G];

        % permute away block structure (x1,x2,x3,y1,y2,y3) -> (x1,y1,x2,y2,x3,y3)
        GGvalues(:,(6*(t-1)+1):(6*t)) = GG(:,[1,4,2,5,3,6]);
      end

      G  = blockfill(2,3,obj.nt, Gvalues);
      GG = blockfill(4,6,obj.nt,GGvalues);

      % Setup element permutation opertor P
      P1 = meshelementpermutation(obj.t,1);
      P2 = meshelementpermutation(obj.t,2);

      obj.GP  =  G*P1;
      obj.GGP = GG*P2;
    end
  end

end

