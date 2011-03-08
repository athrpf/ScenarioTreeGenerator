function A = treeDEconstraint_matrix(ni)
% creates the constraint matrix for the stagewise MILP 

    nnz1 = ni^2;
    nnz2 = 2*ni^2;
    nnz3 = ni;

    nnz = nnz1+nnz2+nnz3;

    irow = zeros(nnz,1);
    jcol = zeros(nnz,1);
    values = zeros(nnz,1);

    n = 1;

    % sum{j in J} eta[i,j] = p[i]
    for ii=1:ni
        for jj=1:ni
            irow(n) = ii;
            jcol(n) = ni*(ii-1)+jj;
            n = n+1;
        end
    end
    values(1:nnz1) = 1;

    % eta[i,j] - z[j] <= 0
    for ii=1:ni
        for jj=1:ni
            irow(n) = ni+ni*(ii-1)+jj;
            jcol(n) = ni*(ii-1)+jj;
            values(n) = 1;
            irow(n+1) = ni+ni*(ii-1)+jj;
            jcol(n+1) = ni^2+jj;
            values(n+1) = -1;
            n = n+2;
        end
    end

    % sum{j in J} z[j] = nc
    for ii=1:ni
        irow(n) = ni+ni^2+1;
        jcol(n) = ni^2+ii;
        n = n+1;
    end
    values(nnz1+nnz2+1:end) = 1;

    A = sparse(irow, jcol, values);
end
