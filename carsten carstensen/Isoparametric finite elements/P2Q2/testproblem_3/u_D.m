function DirichletRandwert=u_Null(u);

% Werte der Funktion auf dem Dirichletrand
% u=[x1,y1; x2,y2; ...]

% der Kondensator hat 3 Dirichlet-Raender: einen aeusseren Rand und 
% 2 Raender an den Kondensator-Platten

% am aeusseren Dirichletrand ist der Wert 0
DirichletRandwert=zeros(size(u,1),1);

% Die Koordinaten relativ zur Mitte des Kondensators
X=u(:,1)-51/2;
Y=u(:,2)-45/2;

% der Dirichlet-Rand an den Kondensator-Platten
InnererRand=find((abs(X)<15).*(abs(Y)<15));

% an der linken Platte ist der Wert 1, an der rechten -1
DirichletRandwert(InnererRand)=sign(X(InnererRand));
