function u = GradLoesung(position)
% Grandient der Loesung des Problems
% position = [x1, y1; x2, y2; ...]

 
[a,r] = cart2pol(position(:,1),position(:,2));
a(find(a<0))=a(find(a<0))+2*pi;
u =  2/3*[-r.^(-1/3).*sin(1/3*a) r.^(-1/3).*cos(1/3*a)];
