%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Mesh3D < Mesh
  %MESH3D

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
  end

  methods(Static)
  end

  methods
    function obj = Mesh3D(p,t,uv)
      %% Constructor
      if nargin < 2,
        [p,t,uv] = readMesh(p);
      end

      obj = obj@Mesh(p,t,uv);

      assert(size(p,1) == 3);

      obj.p   = p;
      obj.p_0 = obj.p;

      obj.nv = size(obj.p,2);
    end

    function p2   = get.p2(obj),     p2 = obj.p(1:2,:);   end
    function p3   = get.p3(obj),     p3 = obj.p;          end
    function p2_0 = get.p2_0(obj), p2_0 = obj.p_0(1:2,:); end
    function p3_0 = get.p3_0(obj), p3_0 = obj.p_0;        end

    function GGP = get.GGP(obj)
      if(isempty(obj.GGP))
        [GP,GGP] = Mesh3D.getGradOps3D(obj);

        obj.GP  = GP;
        obj.GGP = GGP;
      end

      GGP = obj.GGP;
    end

    function GP = get.GP(obj)
      if(isempty(obj.GP))
        [GP,GGP] = Mesh3D.getGradOps3D(obj);

        obj.GP  = GP;
        obj.GGP = GGP;
      end

      GP = obj.GP;
    end

    function GP_3 = get.GP_3(obj)
       GP_3 = obj.GP;
    end

    function h = draw(obj)
      h = trisurf(obj.t', ...
                  obj.p(1,:), ...
                  obj.p(2,:), ...
                  obj.p(3,:), ...
                  zeros(1,obj.nv));
    end
  end

  methods(Static, Access = public)
    function [GP,GGP] = getGradOps3D(mesh)
      %% Compute combined 3D (deformation) gradient operators GP and GGP
      %  (the GGP operator yields 'transposed' deformation gradients)

      % Setup per element deformation gradient operator
      Gvalues  = zeros(3*3,mesh.nt);
      GGvalues = zeros(9*9,mesh.nt);

      n0 = mesh.n0;

      for t=1:mesh.nt,
        % solve 3x3 system for triangle gradient operator
        %   and permute into solution

        p_0 = mesh.p3_0(:,mesh.t(:,t));

        % Permuted Botsch et al. - Deformation
        %   Transfer for Detail-Preserving Surface Editing
        % -> triangle gradient operator

        G = [(p_0(:,2)-p_0(:,1))'; (p_0(:,3)-p_0(:,1))'; n0(:,t)'] \ ...
            [[-1;-1], eye(2); zeros(1,3)];

        Gvalues(:,t) = reshape(G,[],1);

        % deformation gradient operator for blocked coordinates
        GG = [G,zeros(3,3),zeros(3,3);
              zeros(3,3),G,zeros(3,3);
              zeros(3,3),zeros(3,3),G];

        % permute away block structure (x1,x2,x3,y1,y2,y3,z1,z2,z3) -> (x1,y1,z1,x2,y2,z2,x3,y3,z3)
        GGvalues(:,t) = reshape(GG(:,reshape([1:3;4:6;7:9],1,[])),[],1);
      end

      G  = blockfill(3,3,mesh.nt, Gvalues);
      GG = blockfill(9,9,mesh.nt,GGvalues);

      % Setup element permutation opertor P
      P1 = meshelementpermutation(mesh.t,1);
      P3 = meshelementpermutation(mesh.t,3);

      GP  = G*P1;
      GGP = GG*P3;
    end
  end

end
