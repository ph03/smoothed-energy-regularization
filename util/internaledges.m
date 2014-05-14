%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ie, iet, iete, ienv, ieoffset, ieoffsetn, C]=internaledges(t,C)
% Get internal edges of triangulation t.
% Input:
% - t triangulation as obtained by LOADOFFMESH
% - C (optional) connectivity matrix
% Output:
% - ie  internal edges 2 x #IE matrix
% - iet triangle indices of internal edges as 2 x #IE matrix (optional)
% - iete number of edge of internal edges in triangle as 2 x #IE matrix (optional)
% - ienv neighbor vertex to first triangle of internal edge 1 x #IE matrix (optional)
% - ieoffset  index number of second triangle vertices in firsts triangle (1..4, 4 = the neighbor) 3 x #IE matrix (optional)
% - ieoffsetn index number of vertex of first triangle that is adjacent to second triangle (1..3) 1 x #IE matrix (optional)
% - C connectivity matrix, #V x #V sparse (optional)

if nargin<2,
  C=sparse(t([1 2 3],:),t([2 3 1],:),1);
end

[ei,ej]=find(tril(C+C')>1);

ie=[ej,ei]';

if nargout < 2
  return;
end

% Search for two triangles of internal edge using java hashmap

nie=size(ie,2);
nt=size(t,2);

e_to_tri = java.util.Hashtable;

tt=[t;t(1,:)];
for j=1:nt
  for e=1:3
    v=[tt(e,j),tt(e+1,j)];
    v=sort(v);

    edgeid=sprintf('%d,%d',v);
    e_to_tri.put(edgeid,[e_to_tri.get(edgeid),j]);
  end
end

iet=zeros(2,nie);

for e=1:nie
  edgeid=sprintf('%d,%d',ie(:,e));
  iet(:,e)=e_to_tri.get(edgeid);
end

if nargout < 3
  return;
end

iete=zeros(2,nie);

for e=1:nie   %each internal edge
  for k=1:2   %each adjacent triangle
    for i=1:3 %each triangle edge
      v=[tt(i,iet(k,e)),tt(i+1,iet(k,e))];
      v=sort(v)';
      if all(v==ie(:,e))
        iete(k,e)=i;
        break;
      end
    end
  end
end

if nargout < 4
  return;
end

ienv=zeros(1,nie);

for e=1:nie   %each internal edge
  ienv(e) = setdiff(t(:,iet(2,e)),t(:,iet(1,e)));
end

if nargout < 5
  return;
end

ieoffset=zeros(3,nie);

for e=1:nie   %each internal edge
  tleft =t(:,iet(1,e));
  tright=t(:,iet(2,e));

  for i=1:3
    idx = find(tleft==tright(i),1);
    if ~isempty(idx)
      ieoffset(i,e) = idx;
    else
      ieoffset(i,e) = 4;
    end
  end
end

if nargout < 6
  return;
end

ieoffsetn=zeros(1,nie);

for e=1:nie   %each internal edge
  tleft =t(:,iet(1,e));
  tright=t(:,iet(2,e));

  for i=1:3
    idx = find(tright==tleft(i),1);
    if isempty(idx)
      ieoffsetn(e) = i;
      break;
    end
  end
end
end
