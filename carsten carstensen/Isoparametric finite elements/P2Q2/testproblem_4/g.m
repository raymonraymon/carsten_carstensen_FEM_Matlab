function Spannung = g(u)

% Spannung auf dem Neumannrand
% u = [x1,y1; x2,y2; ...]
 
 [alpha,r]=cart2pol(u(:,1),u(:,2));
 alpha(find(alpha<0))=alpha(find(alpha<0))+2*pi;
 Spannung=2/3*r.^(-1/3).*sin(2/3*alpha);
