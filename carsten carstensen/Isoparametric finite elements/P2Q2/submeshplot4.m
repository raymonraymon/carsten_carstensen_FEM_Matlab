function h = submeshplot4(coordinates, elements, u, granularity)
%SUBMESHPLOT4 Plot the given Q2 quadrilaterals with a P1 submesh
%  SUBMESHPLOT4(coordinates, elements, u, granularity)
%  COORDINATES coordinates of points, NPx2
%  ELEMENTS list of points with 4 vertices, 4 edgenodes and one center,
%    NEx9
%  U solution on the points in hierarchical basis
%  GRANULARITY number of subdivisions

 if isempty(elements)
   return
 end
% generate coordinates on the reference quadrilateral
 [Y,X] = meshgrid(-granularity:2:granularity, -granularity:2:granularity);
 sm_coords_ref = [X(:), Y(:)]/granularity;
% generate triangles on the reference quadrilateral
% as patch doesn't interpolate nicely, have P1 on the submesh
 N = granularity + 1;
 pnts = reshape(1:N*N, N, N);
 pnts_ll = pnts(1:end-1, 1:end-1); %% lower left
 pnts_lr = pnts(1:end-1, 2:end); %% lower right
 pnts_ul = pnts(2:end, 1:end-1); %% upper left
 pnts_ur = pnts(2:end, 2:end); %% upper right
 sm_elems = [pnts_ll(:), pnts_ul(:), pnts_lr(:); ...
            pnts_ul(:), pnts_ur(:), pnts_lr(:)];
% generate the patches for each triangle and interpolate solution
 vertices = [];
 coords = [];
 U = [];
 inc = size(sm_coords_ref,1);
 pm = 1 - sm_coords_ref;
 pp = 1 + sm_coords_ref;
 p2 = 1 - sm_coords_ref.^2;
 psi = [pm(:,1).*pm(:,2), pp(:,1).*pm(:,2), ...
        pp(:,1).*pp(:,2), pm(:,1).*pp(:,2)]/4;
 psi = [psi, [p2(:,1).*pm(:,2), p2(:,2).*pp(:,1), ...
              p2(:,1).*pp(:,2), p2(:,2).*pm(:,1)]/2, p2(:,1).*p2(:,2)];
% compute offsets on edges
 edgeOff = zeros(4,2,size(elements,1));
 ind1 = find(elements(:,5:8));
 if ~isempty(ind1)
   [r,c] = ind2sub([size(elements,1),4],ind1);
   ind2 = sub2ind([size(elements,1),4], r, rem(c,4)+1);
   indM = ind1 + 4*size(elements,1);
   tmp = coordinates(elements(indM),:) - ...
         (coordinates(elements(ind1),:) + coordinates(elements(ind2),:))/2;
   ind = sub2ind([4,size(elements,1)*2],c,2*r-1);
   edgeOff(ind) = tmp(:,1);
   edgeOff(ind+4) = tmp(:,2);
 end
% compute offsets on the center
 centerOff = zeros(size(elements,1),2);
 ind = find(elements(:,9));
 if ~isempty(ind)
   contEdge = reshape(sum(edgeOff(:,:,ind),1)/2, 2, length(ind))';
   linearmid = mean(reshape(coordinates(elements(ind,1:4),:), ...
                            [length(ind),4,2]), 2);
   centerOff(ind,:) = coordinates(elements(ind,9),:) - ...
       reshape(linearmid, length(ind), 2) - contEdge;
 end
% assemble the submeshes
 for n = 1:size(elements,1)
   vertices = [vertices; sm_elems + inc*(n-1)];
   K_T = find(elements(n,:));
   uloc = zeros(9,1);
   uloc(K_T) = u(elements(n,K_T));
   coords = [coords; psi*[coordinates(elements(n,1:4),:); ...
                       edgeOff(:,:,n); centerOff(n,:)]];
   U = [U; psi*uloc];
 end
% plot
 col = mean(U(vertices),2);
 hh = trisurf(vertices, coords(:,1), coords(:,2), U, col, 'edgecolor','none');
 if nargout
   h = hh;
 end
