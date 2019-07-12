function h = drawgrid(coordinates, elements3, elements4, u, granularity)
%SUBMESHPLOT4 Draw the mesh for the given P2/Q2 triangulation
%  DRAWGRID(coordinates, elements3, elements4, u, granularity)
%  COORDINATES coordinates of points, NPx2
%  ELEMENTS3 list of points with 3 vertices and 3 edgenodes NEx6
%  ELEMENTS4 list of points with 4 vertices, 4 edgenodes and one center,
%    NEx9
%  U solution on the points in hierarchical basis
%  GRANULARITY number of subdivisions

% create the edges
 if ~isempty(elements3)
   vert = reshape([elements3(:,1:3), elements3(:,[2,3,1])], ...
                  3*size(elements3,1), 2);
   middle = reshape(elements3(:,4:6), 3*size(elements3,1), 1);
 else
   vert = [];
   middle = [];
 end
 if ~isempty(elements4)
   vert = [vert; reshape([elements4(:,1:4), elements4(:,[2,3,4,1])], ...
                  4*size(elements4,1), 2)];
   middle = [middle; reshape(elements4(:,5:8), 4*size(elements4,1), 1)];
 end
 [verts,I] = unique(sort(vert, 2),'rows');
 mids = middle(I);
% curved edges
 I = find(mids);
 if ~isempty(I)
   offset = coordinates(mids(I),:) - ...
            (coordinates(verts(I,1),:) + coordinates(verts(I,2),:))/2;
   l = (0:granularity)/granularity;
   lx = coordinates(verts(I,1),1)*l + ...
        coordinates(verts(I,2),1)*(1-l) + offset(:,1)*((1-l).*l)*4;
   ly = coordinates(verts(I,1),2)*l + ...
        coordinates(verts(I,2),2)*(1-l) + offset(:,2)*((1-l).*l)*4;
   U = u(verts(I,1))*l + u(verts(I,2))*(1-l) + u(mids(I))*((1-l).*l)*4;
   hh = plot3(lx', ly', U', 'k-');
 else
   hh = [];
 end
 hld = ishold;
 hold on;
% linear edges
 I = find(~mids);
 if ~isempty(I)
   lx = reshape(coordinates(verts(I,:),1), length(I), 2);
   ly = reshape(coordinates(verts(I,:),2), length(I), 2);
   U = reshape(u(verts(I,:)), length(I), 2);
   hh = [hh; plot3(lx', ly', U', 'k-')];
 end
 if ~ishold
   hold off;
 end
 if nargout
   h = hh;
 end
