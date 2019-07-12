function DirichletRandwert = u_Null(u)
% Werte der Funktion auf dem Dirchlet-Rand
% u = [x1, y1; x2, y2; ...]
 
 [alpha,r] = cart2pol(u(:,1),u(:,2));
 alpha(find(alpha<0))=alpha(find(alpha<0))+2*pi;
 DirichletRandwert=r.^(2/3).*sin(2/3*alpha);
