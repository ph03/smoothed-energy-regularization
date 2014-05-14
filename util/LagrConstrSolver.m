%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef LagrConstrSolver < handle
  %LAGRCONSTSOLVER

  properties
    m
    n

    L
    U

    p
    q

    b = []
  end

  methods
    function obj = LagrConstrSolver(E,A_or_cidxs,b)
      %% Constructor
      %  Minimizes f(x) = 1/2*x'*E*x - x'*b + d s.t. A*x=c using
      %  Lagrange multipliers.
      %  If cidxs are given A is a corresponding identity matrix.

      if nargin > 2 && ~isempty(b),
        obj.b=b;
      end

      obj.m=size(E,1);
      obj.n=size(E,2);

      if issparse(A_or_cidxs),
        A=A_or_cidxs;
        nconstr = size(A,1);
      else
        assert(size(A_or_cidxs,1)<=1);
        nconstr=size(A_or_cidxs,2);
        A=sparse(1:nconstr,A_or_cidxs,1,nconstr,obj.n);
      end

      EE=[E,A';A,sparse(nconstr,nconstr)];

      [obj.L,obj.U,obj.p,obj.q]=lu(EE,'vector');
    end

    function [x,lambda] = solve(obj,c,b)
      if nargin < 3,
        assert(~isempty(obj.b));
        b=obj.b;
      end

      bb=[b;c];

      x(obj.q,:)=obj.U\(obj.L\bb(obj.p,:));

      if nargout > 1,
        lambda=x((obj.m+1):end,:);
      end

      x((obj.m+1):end,:)=[];
    end

    function setRHS(obj,b)
      obj.b=b;
    end
  end
end
