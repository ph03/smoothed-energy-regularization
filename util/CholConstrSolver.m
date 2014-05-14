%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef CholConstrSolver < handle
  %CHOLCONSTRSOLVER

  properties
    m
    n

    cidxs
    fidxs

    L00
    A01

    b = []

    cholPerm
    cholPermInv
  end

  methods
    function obj = CholConstrSolver(A,cidxs,b)
      %% Constructor
      assert(size(cidxs,1) <= 1);

      nconstr = size(cidxs,2);

      obj.cidxs = cidxs;

      obj.m=size(A,1);
      obj.n=size(A,2);

      obj.fidxs=setdiff(1:obj.n,cidxs);

      A00    =A(obj.fidxs,obj.fidxs);
      obj.A01=A(obj.fidxs,obj.cidxs);

      [obj.L00,bCholPosDef,obj.cholPerm]=chol(A00,'lower','vector');
      obj.cholPermInv(obj.cholPerm)     =1:(obj.n-nconstr);

      assert(bCholPosDef == 0);

      if nargin > 2 && ~isempty(b),
        obj.b=b;
      end
    end

    function x = solve(obj,constr,b)
      if nargin < 3,
        assert(~isempty(obj.b));
        b=obj.b;
      end

      b(obj.cidxs,:)=[];
      b=b-obj.A01*constr;

      y=obj.L00'\(obj.L00\b(obj.cholPerm,:));
      y=y(obj.cholPermInv,:);

      x=zeros(obj.m,size(b,2));
      x([obj.fidxs, obj.cidxs],:) = [y;constr];
    end

    function setRHS(obj,b)
      obj.b=b;
    end
  end

end
