function v=rand_unit_vec(n)
%Generate a random vector from uniform distribution on the unit sphere B(n)
v0=randn(n,1);
v=v0/norm(v0);
end