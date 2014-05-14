%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S=blockfill(bm,bn,blknum,v)
%BLOCKFILL Returns the sparse matrix S which is filled blockwise in
% column major order with the values v. Parameters are the same as for
% blockfillidx.
%
% See blockfillidx.

[i,j]=blockfillidx(bm,bn,blknum);

S=sparse(i,j,v);
end
