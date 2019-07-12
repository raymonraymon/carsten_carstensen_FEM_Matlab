function P2Q2Iso2D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 str = './testproblem_7/';
 addpath(str);
 granularity = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Initialize
 load coordinates.dat;
 eval('load elements3.dat','elements3 = [];'); 
 eval('load elements4.dat','elements4 = [];'); 
 load Dirichlet.dat;
 eval('load Neumann.dat','Neumann = [];'); 
 A = sparse(size(coordinates,1),size(coordinates,1));
 b = zeros(size(coordinates,1),1); u = b; v = b;
% Local stiffness matrix and volume forces for elements with three vertices
 [psi,psi_r,psi_s,kappa] = quad3(7);
 N = [1,1,0;0,1,1;1,0,1]/2;
 for j = 1 : size(elements3,1)
   K_T = find(elements3(j,:));
   P = zeros(6,2);
   P(K_T,:) = coordinates(elements3(j,K_T),:);
   P(4:6,:) = P(4:6,:) + ((elements3(j,4:6)==0)' * [1,1]) .* (N * P(1:3,:));
   D(1:3,:) = P(1:3,:);
   D(4:6,:) = P(4:6,:) - (N * P(1:3,:));
   M = zeros(6,6);
   for m = 1 : size(kappa,2)
     D_Psi = [psi_r(m,:);psi_s(m,:)] * D;
     F = inv(D_Psi) * [psi_r(m,:);psi_s(m,:)];
     det_D_Psi(m) = abs(det(D_Psi));
     M = M + kappa(m) * (F' * F) * det_D_Psi(m);
   end
   A(elements3(j,K_T),elements3(j,K_T)) = ...
       A(elements3(j,K_T),elements3(j,K_T)) + M(K_T,K_T);
   d = kappa .* det_D_Psi .* f(psi * D)' * psi;
   b(elements3(j,K_T)) = b(elements3(j,K_T)) + d(K_T)';
 end
% Local stiffness matrix and volume forces for elements with four vertices
 [phi,phi_xi,phi_eta,gamma] = quad4(9);
 K = [1,1,0,0;0,1,1,0;0,0,1,1;1,0,0,1]/2;        
 L = [-1,-1,-1,-1,2,2,2,2]/4;                    
 for j = 1 : size(elements4,1)
   J_T = find(elements4(j,:));                
   P = zeros(9,2);
   P(J_T,:) = coordinates(elements4(j,J_T),:);    
   P(5:8,:) = P(5:8,:) + ((elements4(j,5:8)==0)' * [1,1]) .* (K * P(1:4,:));  
   P(9,:) = P(9,:) + ((elements4(j,9)==0)' * [1,1]) .* (L * P(1:8,:)); 
   C(1:4,:) = P(1:4,:);
   C(5:8,:) = P(5:8,:) - (K * P(1:4,:));         
   C(9,:) = P(9,:) - (L * P(1:8,:)); 
   M = zeros(9,9);
   for m = 1 : size(gamma,2)
     D_Phi = [phi_xi(m,:);phi_eta(m,:)] * C;
     F = inv(D_Phi) * [phi_xi(m,:);phi_eta(m,:)];
     det_D_Phi(m) = abs(det(D_Phi));
     M = M + gamma(m) * (F' * F) * det_D_Phi(m);
   end 
   A(elements4(j,J_T),elements4(j,J_T)) = ...
       A(elements4(j,J_T),elements4(j,J_T)) + M(J_T,J_T);
   d = gamma .* det_D_Phi .* f(phi * C)' * phi;
   b(elements4(j,J_T)) = b(elements4(j,J_T)) + d(J_T)';
 end
% Neumann conditions
 [phi_E,phi_E_dt,delta_E] = quadN(3);
 for j = 1 : size(Neumann,1)
   J_E = find(Neumann(j,:));
   P = zeros(3,2);    
   P(J_E,:) = coordinates(Neumann(j,J_E),:);
   P(3,:) = P(3,:) + ((Neumann(j,3)==0)' * [1,1]) .* (P(1,:) + P(2,:))/2;
   G(1:2,:) = P(1:2,:);
   G(3,:) = P(3,:) - (P(1,:) + P(2,:))/2;
   norm_Phi_E_dt = sqrt(sum((phi_E_dt * G)'.^2));
   d = delta_E .* g(phi_E * G)' .* norm_Phi_E_dt * phi_E;
   b(Neumann(j,J_E)) = b(Neumann(j,J_E)) + d(J_E)';
 end
% Dirichlet conditions
 ind = find(Dirichlet(:,3));
 u(unique(Dirichlet(:,1:2))) = u_D(coordinates(unique(Dirichlet(:,1:2)),:));
 u(Dirichlet(ind,3)) = u_D(coordinates(Dirichlet(ind,3),:)) - ...
     (u(Dirichlet(ind,1)) + u(Dirichlet(ind,2)))/2;
 b = b - A * u;
% Hanging nodes
 eval('load hanging_nodes.dat;hn=hanging_nodes;','hn = [];'); 
 if ~isempty(hn)
   M = [1,1,2,-2,0,0;1,1,3,-2,-4,0;1,1,3,-2,0,-4];
   B = sparse(3*size(hn,1), size(coordinates,1));
   for j = 1:size(hn,1)
     B((1:3)+(j-1)*3,hn(j,:)) = M;
   end
   lambdas = size(coordinates,1)+(1:3*size(hn,1));
   A = [A,B';B,sparse(3*size(hn,1), 3*size(hn,1))];
   b = [b;zeros(3*size(hn,1),1)];
   v = [v;zeros(3*size(hn,1),1)];
 else
   lambdas = [];
 end
% Compute solution in free nodes
 freeNodes = [setdiff(1:size(coordinates,1),unique(Dirichlet)), lambdas];
 v(freeNodes) = A(freeNodes,freeNodes) \ b(freeNodes);
 if ~isempty(hn)
   v(size(coordinates,1)+1:end,:) = [];
 end
% Display solution
 submeshplot3(coordinates, elements3, v+u, granularity);
 hold on
 submeshplot4(coordinates, elements4, v+u, granularity);
 drawgrid(coordinates, elements3, elements4, v+u, granularity);
 hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 rmpath(str);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [psi,psi_r,psi_s,kappa] = quad3(K_3);
 switch K_3
   case 1
     r = 1/3;
     s = 1/3;
     kappa = 1/2;
   case 3
     r = [1,4,1]'/6;
     s = [1,1,4]'/6;
     kappa = [1,1,1]/6;
   otherwise     
     pos = [6-sqrt(15),9+2*sqrt(15),6+sqrt(15),9-2*sqrt(15),7]/21;
     r = pos([1,2,1,3,3,4,5])';
     s = pos([1,1,2,4,3,3,5])';
     wts = [155-sqrt(15),155+sqrt(15),270]/2400;
     kappa = wts([1,1,1,2,2,2,3]);
 end   
 one = ones(size(kappa,2),1);
 psi = [1-r-s,r,s,4*r.*(1-r-s),4*r.*s,4*s.*(1-r-s)];
 psi_r = [-one,one,0*one,4*(1-2*r-s),4*s,-4*s];
 psi_s = [-one,0*one,one,-4*r,4*r,4*(1-r-2*s)];

function [phi,phi_xi,phi_eta,gamma] = quad4(K_4);
 switch K_4
   case 1 
     xi = 0;
     eta = 0;
     gamma = 4;
   case 4
     xi = sqrt(1/3) * [-1,1,1,-1]';
     eta = sqrt(1/3) * [-1,-1,1,1]';
     gamma = [1,1,1,1];
   otherwise
     xi = sqrt(3/5) * [-1,0,1,-1,0,1,-1,0,1]';   
     eta = sqrt(3/5) * [-1,-1,-1,0,0,0,1,1,1]';
     gamma = [25,40,25,40,64,40,25,40,25]/81;     
 end
 phi = [(1-xi).*(1-eta)/2,(1+xi).*(1-eta)/2,...
       (1+xi).*(1+eta)/2,(1-xi).*(1+eta)/2,...
       (1-xi.^2).*(1-eta),(1+xi).*(1-eta.^2),...
       (1-xi.^2).*(1+eta),(1-xi).*(1-eta.^2),...
       2*(1-xi.^2).*(1-eta.^2)]/2;
 phi_xi = [-(1-eta)/2,(1-eta)/2,(1+eta)/2,-(1+eta)/2,...
       -2*xi.*(1-eta),1-eta.^2,-2*xi.*(1+eta),-1+eta.^2,...
       -4*xi.*(1-eta.^2)]/2; 
 phi_eta = [-(1-xi)/2,-(1+xi)/2,(1+xi)/2,(1-xi)/2,...
       -1+xi.^2,-2*(1+xi).*eta,1-xi.^2,-2*(1-xi).*eta,...
       -4*(1-xi.^2).*eta]/2; 

function [phi_E,phi_E_dt,delta_E] = quadN(K_N);
 switch K_N
   case 1
     t = 0;
     delta_E = 2;     
   otherwise
     t = sqrt(3/5) * [-1,0,1]';
     delta_E = [5,8,5]/9;
 end     
 one = ones(size(delta_E,2),1);
 phi_E = [1-t,1+t,2*(1-t).*(1+t)]/2;
 phi_E_dt = [-one,one,-4*t]/2;
