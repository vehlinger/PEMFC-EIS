function C = freq_response(C,n,nj,params,op_cond,C_ss)

    jcount = 0;          % current iteration
    dC = 1e-8*ones(1,n); % Delta C = small variation in value of C

    rtol = 1e-6;
    atol = 1e-9;
    kerr = 1;
    kerrg = 1;

    itmax = 10;
    
    while (kerr == 1 || kerrg == 1) && jcount < itmax
        jcount = jcount+1;   % update iteration
        CC = C;              % initialize CC

        C = autoband_f(n,nj,C,dC,params,op_cond,C_ss);

        kerr = 0;
        kerrg = 0;
        for j = 1:nj
            for i = 1:n
                if kerr == 0 && kerrg == 0
                    if abs(C(j,i)) > rtol*abs(CC(j,i))
                        kerr = 1;
                    end
                    if kerr == 1 && abs(abs(C(j,i))<atol)
                        kerr = 0;
                    end
                end
            end
            for i = 1:n
                C(j,i) = CC(j,i)+C(j,i);
            end
        end
    end
end
