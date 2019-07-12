function plotisoparametric(c4n,n4e,u)
if nargin < 3 | isempty(u), u = zeros(size(c4n,1),1); end;
hold on;
for j = 1:size(n4e,1)
 element = [];
 for k = 1:fix(size(n4e,2)/2)
  y = [c4n(n4e(j,k),:),u(n4e(j,k));
       c4n(n4e(j,mod(k,fix(size(n4e,2)/2))+1),:),...
       u(n4e(j,mod(k,fix(size(n4e,2)/2))+1))];
  y(3,:) = eval(strcat('[c4n(n4e(j,k+fix(size(n4e,2)/2)),:)',...
           ',u(n4e(j,k+fix(size(n4e,2)/2)))+mean(y(1:2,3))]'),'mean(y(1:2,:))');
  t = linspace(-1,1,100);
  element = [element,...
             y(1,:)'*(1-t)/2+y(2,:)'*(1+t)/2+(y(3,:)'-mean(y(1:2,:))')*((1-t).*(1+t))];
 end
 plot3(element(1,:),element(2,:),element(3,:),'k','EraseMode','background');
end
hold off; 
