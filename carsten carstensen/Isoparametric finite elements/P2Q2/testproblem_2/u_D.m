function rueck = u_Null(x);

rueck = ones(size(x,1),1);
rueck(find(x(:,1)>27.5)) = -ones(size(find(x(:,1)>27.5)));
