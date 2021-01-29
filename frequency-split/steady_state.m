function C = steady_state(C,n,nj,params,op_cond)
    jcount = 0;          % current iteration
    dC = 1e-8*[ones(1,4) 1e-4 1]; % Delta C = small variation in value of C

    rtol = 1e-6;
    atol = 1e-9;
    kerr = 1;
    kerrg = 1;

    itmax = 10;
    
    while (kerr == 1 || kerrg == 1) && jcount < itmax
        jcount = jcount+1;   % update iteration
        CC = C;              % initialize CC

        C = autoband(n,nj,C,dC,params,op_cond,[],0);

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
    
    L = op_cond(1);
    x = (0:L/(nj-1):L)*10000;

%     if plot_on == 1
%         figure(1)
%         subplot(1,3,1)
%         plot(x,C(:,1))
%         hold on
%         plot(x,C(:,3))
%         hold off
%         xlabel('L (\mum)')
%         ylabel('Current Density (A/cm^2)')
%         tstr = sprintf('t = 0 s');
%         title(tstr)
%         legend('i_1','i_2')
%         subplot(1,3,2)
%         plot(x,C(:,2))
%         hold on
%         plot(x,C(:,4))
%         hold off
%         legend('v_1','v_2')
%         xlabel('L (\mum)')
%         ylabel('Potential (V)')
%         subplot(1,3,3)
%         yyaxis left
%         plot(x,C(:,5))
%         ylabel('Concentration (mol/cm^3)')
%         yyaxis right
%         plot(x,C(:,6))
%         ylabel('Flux (mol/cm^2/s)')
%         xlabel('L (\mum)')
%     end

end
