function eq = eqn_fsplit(j,jp,k,dC,C,nj,params,op_cond,C_ss)
    C(jp,k) = C(jp,k)+dC;
    % unknowns at each point
    ii1 = 1; iv1 = 2; ii2 = 3; iv2 = 4; iyO2 = 5; iNO2 = 6;
    n = params(1);
    % parameters
    sigma = params(4); kappa = params(5); D = params(8);
    % operating conditions
    L = op_cond(1);   deltaV = op_cond(3);  dx = L/(nj-1);
    
    %% Equation 1: Charge Balance
    % Real
    if j == 1
        eq(ii1) = C(j,ii1);
    else
        eq(ii1) = (C(j,ii1)-C(j-1,ii1))/dx+(C(j,ii2)-C(j-1,ii2))/dx;
    end
    % Imaginary
    if j == 1
        eq(ii1+n) = C(j,ii1+n);
    else
        eq(ii1+n) = (C(j,ii1+n)-C(j-1,ii1+n))/dx+...
            (C(j,ii2+n)-C(j-1,ii2+n))/dx;
    end
    
    %% Equation 2: Ohm's Law
    % Real
    if j == nj 
        eq(iv1) = C(j,iv1) - deltaV;
    else
        eq(iv1) = C(j,ii1) + sigma*(C(j+1,iv1)-C(j,iv1))/dx;
    end
    % Imaginary
    if j == nj
        eq(iv1+n) = C(j,iv1+n);
    else
        eq(iv1+n) = C(j,ii1+n) +...
            sigma*(C(j+1,iv1+n)-C(j,iv1+n))/dx;
    end
    
    %% Equation 3: Flux (no diffusion or convection)
    % Real
    if j < nj
        eq(ii2) = C(j,ii2) + kappa*(C(j+1,iv2)-C(j,iv2))/dx;
    else
        eq(ii2) = C(j,ii2);
    end
    % Imaginary 
    if j < nj
        eq(ii2+n) = C(j,ii2+n) +...
            kappa*(C(j+1,iv2+n)-C(j,iv2+n))/dx;
    else
        eq(ii2+n) = C(j,ii2+n);
    end
    
    %% Equation 4: Polarization (kinetics)
    % Real
    if j == 1
        eq(iv2) = C(j,iv2);
    else
        eq(iv2) = fluxleft_fsplit(j,ii2,C,C_ss,params,op_cond,dx)-...
            fluxright_fsplit(j,ii2,C,C_ss,params,op_cond,dx);
    end
    % Imaginary
    if j == 1
        eq(iv2+n) = C(j,iv2+n);
    else
        eq(iv2+n) = fluxleft_fsplit(j,ii2+n,C,C_ss,params,op_cond,dx)-...
            fluxright_fsplit(j,ii2+n,C,C_ss,params,op_cond,dx);
    end
    
    %% Equation 5: Concentration Gradient (Fick's Law)
    % Real
    if j == 1    
        eq(iyO2) = fluxright_fsplit(j,iNO2,C,C_ss,params,op_cond,dx);
    else
        R = 83.14; % cm3 bar / mol K
        T = op_cond(2); p = op_cond(5);
        CT = p/(T*R);
        eq(iyO2) = C(j,iNO2) + D*CT*(C(j,iyO2)-C(j-1,iyO2))/dx;
    end
    % Imaginary
    if j == 1    
        eq(iyO2+n) = fluxright_fsplit(j,iNO2+n,C,C_ss,params,op_cond,dx);
    else
        R = 83.14; % cm3 bar / mol K
        T = op_cond(2); p = op_cond(5);
        CT = p/(T*R);
        eq(iyO2+n) = C(j,iNO2+n) + D*CT*(C(j,iyO2+n)-C(j-1,iyO2+n))/dx;
    end    
    %% Equation 6: Flux (conservation of mass)
    % Real
    if j < nj
        eq(iNO2) = fluxleft_fsplit(j,iNO2,C,C_ss,params,op_cond,dx)-...
            fluxright_fsplit(j,iNO2,C,C_ss,params,op_cond,dx);
    else
        eq(iNO2) = C(j,iyO2); 
    end
    % Imaginary
    if j < nj
        eq(iNO2+n) = fluxleft_fsplit(j,iNO2+n,C,C_ss,params,op_cond,dx)-...
            fluxright_split(j,iNO2+n,C,C_ss,params,op_cond,dx);
    else
        eq(iNO2+n) = C(j,iyO2+n); 
    end
end