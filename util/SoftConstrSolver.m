%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef SoftConstrSolver < handle
  %SOFTCONSTRSOLVER

  properties
    scidxs
    hcidxs

    scindicators
    hcindicators

    b = []

    w
    C

    csolver
  end

  methods
    function obj = SoftConstrSolver(A,soft_cidxs_or_soft_hard_cidxs,b,w,csolverf)
      %% Constructor
      % This is an adaptor for the other (hard constraint) ConstrSolvers
      % to support soft and mixed soft / hard constraints.
      %
      % If soft_cidxs_or_soft_hard_cidxs is a vector the indices are
      % used as 'soft constraints' by regularization with weight w (default 1).
      %
      % If soft_cidxs_or_soft_hard_cidxs is a cell array the first indices are
      % used as 'soft constraints' by regularization with weight w (default
      % 1) and the second indices are used as explicit 'hard constraints'
      %
      % A ConstrSolvers can be chosen by csolverf
      % (defaults to @(A,cidxs,b)CholConstrSolver(A,cidxs,b))

      if iscell(soft_cidxs_or_soft_hard_cidxs),
        %% Mixed soft / hard constraints
        assert(numel(soft_cidxs_or_soft_hard_cidxs) == 2);

        obj.scidxs = soft_cidxs_or_soft_hard_cidxs{1};
        obj.hcidxs = soft_cidxs_or_soft_hard_cidxs{2};

        assert(size(obj.hcidxs,1) == 1);
      else
        %% Soft constraints only
        obj.scidxs = soft_cidxs_or_soft_hard_cidxs;
        obj.hcidxs = [];
      end

      assert(size(obj.scidxs,1) == 1);
      assert(isempty(intersect(obj.scidxs,obj.hcidxs)));

      if nargin > 2 && ~isempty(b),
        obj.b=b;
      end

      if nargin < 4,
        obj.w = 1;
      else
        obj.w = w;
      end

      if nargin < 5,
        csolverf = @(A,cidxs,b)CholConstrSolver(A,cidxs,b);
      end

      nsconstr = size(obj.scidxs,2);
      nhconstr = size(obj.hcidxs,2);

      %% soft constraints selector matrix C, adapt / regularize system
      obj.C = sparse(1:nsconstr,obj.scidxs,1,nsconstr,size(A,2));
      A     = A + obj.w * obj.C'*obj.C;

      %% use a ConstrSolver for actual solving with hard constraints
      obj.csolver = csolverf(A,obj.hcidxs,obj.b);

      %% setup constraint indicators
      [~,isort]        =  sort([obj.scidxs,obj.hcidxs]);
      obj.hcindicators =  [false(1,nsconstr),true(1,nhconstr)];
      obj.hcindicators =  obj.hcindicators(isort);
      obj.scindicators = ~obj.hcindicators;
    end

    function x = solve(obj,constr,b)
      if nargin < 3,
        assert(~isempty(obj.b));
        b=obj.b;
      end

      %% adapt rhs with soft constraint terms
      b=b+obj.w*(obj.C'*constr(obj.scindicators,:));

      %% solve with hard constraints
      x=obj.csolver.solve(constr(obj.hcindicators,:),b);
    end

    function setRHS(obj,b)
      obj.b=b;
    end
  end
end
