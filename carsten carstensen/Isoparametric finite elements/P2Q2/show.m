function show(coordinates, elements3, elements4, x);
% Process elements with four vertices 
K = [1,1,0,0;0,1,1,0;0,0,1,1;1,0,0,1]/2;        
L = [-1,-1,-1,-1,2,2,2,2]/4;                    
for j = 1 : size(elements4,1)
  J_T = find(elements4(j,:)); 
  P = zeros(9,2);
  P(J_T,:) = coordinates(elements4(j,J_T),:);    
  P(5:8,:) = P(5:8,:) + ((elements4(j,5:8)==0)' * [1,1]) .* (K * P(1:4,:));  
  P(9,:) = P(9,:) + ((elements4(j,9)==0)' * [1,1]) .* (L * P(1:8,:)); 
  n_ind = find(elements4(j,:)==0); 
  elements4(j,n_ind) = size(coordinates,1) + [1:size(n_ind,2)];
  coordinates = [coordinates;P(n_ind,:)];
end 
x = [x;zeros(size(coordinates,1)-size(x,1),1)];
fractions = [4,0,0,0,0,0,0,0,0;0,4,0,0,0,0,0,0,0;...
             0,0,4,0,0,0,0,0,0;0,0,0,4,0,0,0,0,0;...
             2,2,0,0,4,0,0,0,0;0,2,2,0,0,4,0,0,0;...
             0,0,2,2,0,0,4,0,0;2,0,0,2,0,0,0,4,0;...
             1,1,1,1,2,2,2,2,4]/4;
for j = 1 : size(elements4,1)
  r(elements4(j,:)) = fractions * x(elements4(j,:));
end
if ~isempty(elements4)
  triangles = [elements4(:,[1,5,9]);elements4(:,[1,9,8]);...
               elements4(:,[5,2,6]);elements4(:,[5,6,9]);... 
               elements4(:,[9,6,3]);elements4(:,[9,3,7]);... 
               elements4(:,[8,9,7]);elements4(:,[8,7,4])]';
end
% Process elements with three vertices 
N = [1,1,0;0,1,1;1,0,1]/2;
for j = 1 : size(elements3,1)
  K_T = find(elements3(j,:));
  P = zeros(6,2);
  P(K_T,:) = coordinates(elements3(j,K_T),:);
  P(4:6,:) = P(4:6,:) + ((elements3(j,4:6)==0)' * [1,1]) .* (N * P(1:3,:));
  n_ind = find(elements3(j,:)==0); 
  elements3(j,n_ind) = size(coordinates,1) + [1:size(n_ind,2)];
  coordinates = [coordinates;P(n_ind,:)];
end 
x = [x;zeros(size(coordinates,1)-size(x,1),1)];
fractions = [2,0,0,0,0,0;0,2,0,0,0,0;0,0,2,0,0,0;...
             1,1,0,2,0,0;0,1,1,0,2,0;1,0,1,0,0,2]/2;
for j = 1 : size(elements3,1)
  r(elements3(j,:)) = fractions * x(elements3(j,:));
end
if ~isempty(elements3)
  triangles = [triangles,[elements3(:,[1,4,6]);elements3(:,[2,5,4]);... 
                      elements3(:,[3,6,5]);elements3(:,[4,5,6])]'];
end
% Display numerical solution
trisurf(triangles',coordinates(:,1),coordinates(:,2),r');
