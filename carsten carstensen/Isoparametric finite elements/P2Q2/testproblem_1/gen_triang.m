function gen_triang(N)

for j = 0 : N
  for k = 0 : N
    coordinates(j*(N+1)+k+1,:) = [j*(N+1)+k+1,k/N,j/N];
  end
end
for j = 1 : N
  for k = 1 : N
    elements((j-1)*N+k,:) = ...
	[(j-1)*N+k,(j-1)*(N+1)+k+[0,1,N+2,N+1],0,0,0,0,0];
  end
end
for j = 1 : N 
  Dirichlet(j,:) = ...
      [j,(N+1)^2-j*(N+1)+1,(N+1)^2-(j+1)*(N+1)+1,0];
end
for j = 1 : N 
  Dirichlet(N+j,:) = [N+j,j,j+1,0];
end
for j = 1 : N 
  Dirichlet(2*N+j,:) = [2*N+j,j*(N+1),(j+1)*(N+1),0];
end
for j = 1 : N 
  Neumann(j,:) = [j,(N+1)^2-[(j-1),j],0];
end

save coordinates.dat coordinates -ASCII
save elements.dat elements -ASCII
save Dirichlet.dat Dirichlet -ASCII
save Neumann.dat Neumann -ASCII
