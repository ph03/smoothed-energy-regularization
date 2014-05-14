%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  ACM TOG - J. Martinez Esturo, C. Rössl, and H. Theisel
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S=blockfill(bm,bn,blknum,v)
%BLOCKFILL Returns the sparse matrix S which is filled blockwise in
% column major order with the values v. Parameters are the same as for
% blockfillidx.
%
% See blockfillidx.

[i,j]=blockfillidx(bm,bn,blknum);

S=sparse(i,j,v);
end
