function plotisoparametric(c4n,n4e)
hold on;
for j = 1:size(n4e,1)
 for k = 1:fix(size(n4e,2)/2)
  y = [c4n(n4e(j,k),:);
       eval('c4n(n4e(j,k+fix(size(n4e,2)/2)),:)','[]');
       c4n(n4e(j,mod(k,fix(size(n4e,2)/2))+1),:)]';
  yy = ppval(spline(1:size(y,2),y),linspace(1,size(y,2),100));
  plot(yy(1,:),yy(2,:),'k');
 end
end
hold off;
