%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = blockreshape(A,blksz,nr,nc)
%% Resize a column vector of quadratic block matrices to a given number
%  of rows and columns.

assert(size(A,2) == blksz);
assert(mod(size(A,1),blksz) == 0);

assert(size(A,1)*size(A,2) == (blksz * nr) * (blksz * nc));

%Engine one line version:
A = reshape(permute(reshape(permute(reshape(A.',blksz,blksz,[]),[2 1 3]), ...
            blksz,blksz*nc,[]),[2 1 3]),blksz*nc,blksz*nr)';

% %Engine in pieces:
% B = reshape(A.',blksz,blksz,[]); %reshape the transpose to 3d with each block a 3×3
% C = permute(B,[2 1 3]); %transpose each block
% D = reshape(C,blksz,blksz*nc,[]); %reshape to nr planes that need to be stacked vertically
% E = permute(D,[2 1 3]); %transpose each plane
% F = reshape(E,blksz*nc,blksz*nr)'; %reshape into a matrix and transpose
end
