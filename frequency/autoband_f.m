function C = autoband_f(n,nj,C,dC,params,op_cond,C_ss)

J = zeros(n*nj);
b = zeros(n*nj,1);

for j = 1:nj
    A = zeros(n,n);     % matrix of dG/dC at j-1 
    B = zeros(n,n);     % matrix of dG/dC at j
    D = zeros(n,n); % matrix of dG/dC at j+1

    % initialize G (k = 1, dC = 0)
    G = eqn_f(j,j,1,0,C,nj,params,op_cond,C_ss); 

    % generate A,B,D matrices
    for k = 1:n
        eq = eqn_f(j,j,k,dC(k),C,nj,params,op_cond,C_ss);
        B(:,k) = -(eq-G)./dC(k);
        if j > 1
            eq = eqn_f(j,j-1,k,dC(k),C,nj,params,op_cond,C_ss);
            A(:,k) = -(eq-G)./dC(k);
        end
        if j < nj
            eq = eqn_f(j,j+1,k,dC(k),C,nj,params,op_cond,C_ss);
            D(:,k) = -(eq-G)./dC(k);
        end
        % construct tridiagonal matrix
        for m = 1:n
            J((m-1)*nj+j,(k-1)*nj+j) = B(m,k);
            if j > 1
                J((m-1)*nj+j,(k-1)*nj+j-1) = A(m,k);
            end
            if j < nj
                J((m-1)*nj+j,(k-1)*nj+j+1) = D(m,k);
            end
        end
        % construct solution vector
        b((k-1)*nj+j) = G(k);
    end
end

Js = sparse(J);
U = Js\b;

C = reshape(U,nj,n);

end