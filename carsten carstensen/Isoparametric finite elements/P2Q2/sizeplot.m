function sizeplot()

f = gcf;
a = gca;
view(2)
yl=get(gca,'ylim');xl=get(gca,'xlim');
axis equal
set(gca,'xlim',xl);set(gca,'ylim',yl);

pbr = get(gca, 'PlotBoxAspectRatio');

r = pbr(1)/pbr(2);
if r < 1
  x = 12;
  y = x/r;
else
  y = 12;
  x = y*r;
end
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 x y]);
p=findobj(gcf,'type','patch');
for i = p'
  v = get(i,'Vertices');
  v(:,3) = 0;
  set(i,'Vertices',v);
end
l=findobj(gcf,'type','line');
for i = l'
  set(i, 'ZData', ones(size(get(i,'ZData'))));
end
