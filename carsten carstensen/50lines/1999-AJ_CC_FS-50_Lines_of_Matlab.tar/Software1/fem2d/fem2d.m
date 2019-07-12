%FEM2D   two-dimensional finite element method for Laplacian.
%
%    FEM2D solves Laplace's equation 
%      - div(grad(u)) = f in Omega
%                   u = u_D on the Dirichlet boundary
%              d/dn u = g   on the Neumann boundary
%    on a geometry described by triangles and parallelograms, respectively
%    and presents the solution graphically.
% 
%    Therefore, FEM2D assembles plane Courant finite elements and calculates
%    a discrete right-hand side as the coefficient vector of the affine
%    element approximation. Volume force and boundary data are given as
%    M-files <f.m>, <g.m>, and <u_d>. FEM2D uses the reduced linear system 
%    of equations to calculate a discrete solution of the Laplace-problem. 
%    The resulting piecewise affine approximation will be graphically 
%    represented.
% 
%    FEM2D loads the mesh data from data-files. The program reads the
%    triangular and/or quadrilateral elements from the files <elements3.dat>
%    and/or <elements4.dat>, respectively. Depending on the mesh one of the
%    two files, but not both, can be omitted. The first column in
%    <elements3.dat> and <elements4.dat> gives the number of each
%    element. This is used for clearness and is not neccesary for the
%    numerical algorithm. The following columns give the number of each
%    node. Nodes of elements are counted anti-clockwise.
% 
%    To adapt the program to a given Laplace equation the user has to
%    specify the data-files <coordinates.dat>, <elements3.dat> and/or
%    <elements4.dat>, <dirichlet.dat>, and <neumann.dat> (optional) and the
%    M-files <f.m>, <u_d.m>, and <g.m> (optional). They have to be in the
%    same directory as <fem2d.m>.
%
%    Remark: This program is a supplement to the paper "Remarks around  
%    50 lines of Matlab: Short finite element implementation" by  
%    J. Alberty, C. Carstensen and S. A. Funken. The reader should 
%    consult that paper for more information.   
%
%
%    M-files you need to run FEM2D
%       <stima3.m>, <stima4.m>, <f.m>, <u_d.m>, <show.m> and <g.m> (optional)
%
%    Data-files you need to run FEM2D
%       <coordinates.dat>, <elements3.dat> and/or <elements4.dat>,
%       <dirichlet.dat>, and <neumann.dat> (optional)

%    J. Alberty, C. Carstensen and S. A. Funken  02-11-99
%    File <fem2d.m> in $(HOME)/acf/fem2d/
%    This program and corresponding data-files give Fig. 3 in 
%    "Remarks around 50 lines of Matlab: Short finite element 
%    implementation"

% Initialisation
clc;clear all;close all;
load coordinates.dat; coordinates(:,1)=[];
eval('load elements3.dat; elements3(:,1)=[];','elements3=[];');
eval('load elements4.dat; elements4(:,1)=[];','elements4=[];');
eval('load neumann.dat; neumann(:,1) = [];','neumann=[];');
load dirichlet.dat; dirichlet(:,1) = [];
FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));
A = sparse(size(coordinates,1),size(coordinates,1));
b = sparse(size(coordinates,1),1);

% Assembly
for j = 1:size(elements3,1)
  A(elements3(j,:),elements3(j,:)) = A(elements3(j,:),elements3(j,:)) ...
      + stima3(coordinates(elements3(j,:),:));
end
for j = 1:size(elements4,1)
  A(elements4(j,:),elements4(j,:)) = A(elements4(j,:),elements4(j,:)) ...
      + stima4(coordinates(elements4(j,:),:));
end

% Volume Forces
for j = 1:size(elements3,1)
  b(elements3(j,:)) = b(elements3(j,:)) + ...
      det([1,1,1; coordinates(elements3(j,:),:)']) * ...
      f(sum(coordinates(elements3(j,:),:))/3)/6;
end
for j = 1:size(elements4,1)
  b(elements4(j,:)) = b(elements4(j,:)) + ...
      det([1,1,1; coordinates(elements4(j,1:3),:)']) * ...
      f(sum(coordinates(elements4(j,:),:))/4)/4;
end

% Neumann conditions
for j = 1 : size(neumann,1)
  b(neumann(j,:))=b(neumann(j,:)) + norm(coordinates(neumann(j,1),:)- ...
      coordinates(neumann(j,2),:)) * g(sum(coordinates(neumann(j,:),:))/2)/2;
end

% Dirichlet conditions 
u = sparse(size(coordinates,1),1);
u(unique(dirichlet)) = u_d(coordinates(unique(dirichlet),:));
b = b - A * u;

% Computation of the solution
u(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);

% graphic representation
show(elements3,elements4,coordinates,full(u));
