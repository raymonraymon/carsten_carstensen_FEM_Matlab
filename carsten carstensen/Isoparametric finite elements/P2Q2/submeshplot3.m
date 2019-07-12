function h = submeshplot3(coordinates, elements, u, granularity)
%SUBMESHPLOT3 Plot the given P2 triangles
%  SUBMESHPLOT3(coordinates, elements, u, granularity)
%  COORDINATES coordinates of points, NPx2
%  ELEMENTS list of points with 3 vertices and 3 edgenodes NEx6
%  U solution on the points in hierarchical basis
%  GRANULARITY number of subdivisions

if isempty(elements)
  return
end

%% generate coordinates on the reference triangle
sm_coords_ref = [];
for n = granularity:-1:0
  sm_coords_ref = [sm_coords_ref, [0:n; (granularity-n)*ones(1, n+1)]];
end
sm_coords_ref = sm_coords_ref'/granularity;

%% generate triangles on the reference triangle
%% as patch doesn't interpolate nicely, have P1 on the submesh
sm_elems = [];
offset = 0;
for n = granularity:-1:1
  bar = [0:n-1; 1:n; (0:n-1)+n];
  lelems = [1:n; 2:n+1; n+2:2*n+1] + offset;
  uelems = [lelems(3,1:end-1);lelems(1, 2:end);lelems(3,2:end)];
  sm_elems = [sm_elems, lelems, uelems];
  offset = offset + n+1;
end

%% generate the patches for each triangle and interpolate solution
vertices = [];
coords = [];
U = [];
inc = size(sm_coords_ref,1);
psi = [1-sum(sm_coords_ref,2), sm_coords_ref(:,1), sm_coords_ref(:,2)];
psi = [psi, ...
       [psi(:,1).*psi(:,2), psi(:,2).*psi(:,3), psi(:,3).*psi(:,1)]*4];
%% compute offsets on edges
centerOff = zeros(3,2,size(elements,1));
ind1 = find(elements(:,4:6));
if ~isempty(ind1)
  [r,c] = ind2sub([size(elements,1),3],ind1);
  ind2 = sub2ind([size(elements,1),3], r, rem(c,3)+1);
  indM = ind1 + 3*size(elements,1);
  tmp = coordinates(elements(indM),:) - ...
        (coordinates(elements(ind1),:) + coordinates(elements(ind2),:))/2;
  ind = sub2ind([3,size(elements,1)*2],c,2*r-1);
  centerOff(ind) = tmp(:,1);
  centerOff(ind+3) = tmp(:,2);
end

for n = 1:size(elements,1)
  vertices = [vertices, sm_elems + inc*(n-1)];
  K_T = find(elements(n,:));
  uloc = zeros(6,1);
  uloc(K_T) = u(elements(n,K_T));
  coords = [coords; psi*[coordinates(elements(n,1:3),:); centerOff(:,:,n)]];
  U = [U; psi*uloc];
end

col = mean(U(vertices'),2);
hh = trisurf(vertices', coords(:,1), coords(:,2), U, col, 'edgecolor','none');
if nargout
  h = hh;
end
