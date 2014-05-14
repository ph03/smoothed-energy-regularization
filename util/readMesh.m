function [p,t,uv] = readMesh(filename)
  % read in vertices and faces from a .off or .obj file
  % Input:
  %   filename  file holding mesh
  % Output:
  %   p  (vertex list)
  %   t  (face list)
  %   uv (texture coordinates list)
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: readOBJ, readOBJfast, readOFF
  %
  uv = [];

  if ~isempty(regexp(filename,'\.off$'))
    [V,F] = readOFF(filename);
  elseif ~isempty(regexp(filename,'\.obj$'))
    if nargout > 2,
      [V,F,UV] = readOBJ(filename);
      uv = UV';
    else
      try
        [V,F] = readOBJfast(filename);
      catch exception
        fprintf('Fast reader failed, retrying with more robust, slower reader\n');
        [V,F] = readOBJ(filename);
      end
    end
  else
    error('Input file must be .off or .obj file.');
  end

  % adapt to our convention
  p = V';
  t = F';
end
