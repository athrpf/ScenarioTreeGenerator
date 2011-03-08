function Dk = kantorovich(xi, p, nu, q, distance_fun)
% Compute the Kantorovich Distance of two sets of scenarios (p,xi), (q,nu)
% USAGE
% Dk = kantorovich(xi, p, nu, q, distance_fun)
% xi, nu : n_scen x n_stages

assert(size(xi,2) == size(nu,2));
assert(length(p)==size(xi,1));
assert(length(q)==size(nu,1));
c = [];
if nargin > 4
    c = pdist2(xi,nu, distance_fun);
else
    c = pdist2(xi,nu);
end

% construct constraint matrix
[ni,nj] = size(c);
nvar = ni*nj;
nnz = nvar*2;
irow = zeros(nnz,1);
jcol = zeros(nnz,1);
values = ones(nnz,1);

nnzidx = 1;
rowidx = 1;
for jj=1:nj
    for ii=1:ni
        irow(nnzidx) = rowidx;
        jcol(nnzidx) = ni*(jj-1)+ii;
        nnzidx = nnzidx + 1;
    end
    rowidx = rowidx+1;
end

for ii=1:ni
    for jj=1:nj
        irow(nnzidx) = rowidx;
        jcol(nnzidx) = ni*(jj-1)+ii;
        nnzidx = nnzidx + 1;
    end
    rowidx = rowidx+1;
end

A = sparse(irow, jcol, values);
b = zeros(ni+nj,1);
b(1:nj) = q(:);
b(nj+1:end) = p(:);
assert(size(A,1)==length(b));

[x, z, status] = clp([], c(:), [], [], A, b, zeros(nvar,1), ...
                     ones(nvar,1));

if status~=0
    disp(['CLP failed to solve Kantorovich min-cost-flow ' ...
          'problem!']);
    Dk = inf;
end
Dk = x'*c(:);
end