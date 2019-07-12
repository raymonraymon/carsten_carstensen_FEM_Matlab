function rueck = u_D(x);
 val = zeros(size(x,1),1);
 ind = find(x(:,1).^2 + x(:,2).^2 < 26^2 + 8);
 val(ind) = sign(x(ind,1)) .* (26 - x(ind,2)) / 52;
 ind = find((x(:,1)-37).^2 + (x(:,2)-63).^2 < 18^2 + 1);
 val(ind) = 1;
 rueck = val;