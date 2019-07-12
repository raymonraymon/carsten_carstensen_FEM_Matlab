function u = Loesung(position)
% Loesung des Problems
% position = [x1, y1; x2, y2; ...]

 
[a,r] = cart2pol(position(:,1),position(:,2));
a(find(a<0))=a(find(a<0))+2*pi;
u =  r.^(2/3).*sin(2/3*a);
