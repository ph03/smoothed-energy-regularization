%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [triNs,triAs]=trinormals(p,t)
% Get normalized triangle normals and triangle areas.

v1=p(:,t(1,:));
v2=p(:,t(2,:));
v3=p(:,t(3,:));

triNs=cross(v2-v1,v3-v1);
triAs=sqrt(sum(triNs.*triNs,1));
triNs=triNs./repmat(triAs,3,1);

if nargout > 1
  triAs=triAs/2;
end

end
