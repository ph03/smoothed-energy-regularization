%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = spblkdiag(cellspdiag)
%SPBLKDIAG creates a sparse matrix A constructed blockwise of sparse
%          matrices in the cellarray cellspdiag

  %% get total nnz
  nnz_total = 0;
  for s=1:length(cellspdiag),
    nnz_total = nnz_total + nnz(cellspdiag{s});
  end

  %% setup coordinates
  Ai = zeros(1,nnz_total);
  Aj = zeros(1,nnz_total);
  Av = zeros(1,nnz_total);
  krows = 0; kcolumns = 0; knnz = 0;
  for s=1:length(cellspdiag),
    [Si,Sj,Sv] = find(cellspdiag{s});
    nnz_S = nnz(cellspdiag{s});

    Ai((1+knnz):(nnz_S+knnz)) = Si + krows;
    Aj((1+knnz):(nnz_S+knnz)) = Sj + kcolumns;
    Av((1+knnz):(nnz_S+knnz)) = Sv;

    krows    = krows + size(cellspdiag{s},1);
    kcolumns = kcolumns + size(cellspdiag{s},2);
    knnz     = knnz + nnz_S;
  end

  %% construct final sparse block diagonal matrix
  A = sparse(Ai,Aj,Av,krows,kcolumns);
end
