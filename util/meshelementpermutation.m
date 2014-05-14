%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = meshelementpermutation(t,n,coefficientwise)
%% MESHELEMENTPERMUTATION Computes sparse operator P (3*nt*n x nv*n) such
%  that n consecutive values per vertex are mapped to 3*n consecutive
%  output values per triangle.
%  If *coefficientwise* = true (default is false) then the outpur per
%  triangle is given as consecutive coefficientwise vectors, else they are
%  given as the original vectors.

  if nargin < 3
    coefficientwise = false;
  end

  assert(size(t,1) == 3);

  nt = size(t,2);

  pj = repmat((reshape(t,1,3,[])-1).*n + 1,[n,1,1]) + ...
               repmat((0:(n-1))',[1,3,nt]);

  if coefficientwise,
    pj = permute(pj, [2 1 3]);
  end

  pj = reshape(pj,[],1);

  P=sparse(1:(3*nt*n),pj,1);
end
