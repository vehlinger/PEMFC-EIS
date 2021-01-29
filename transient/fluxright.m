function NR = fluxright(j,i,C,Cp,params,op_cond,dt,dx)
% Calculates the flux exiting the box to the right of point j
% Variable Identifiers
    ii1 = 1; iv1 = 2; ii2 = 3; iv2 = 4; iyO2 = 5; iNO2 = 6;
    n = params(1); nj = params(2);
        
% Flux in the box to the right
if i == ii2
    flux = C(j,i);
elseif i == iNO2
    flux = C(j+1,i);
end

% Reaction terms
F = 96845; % C/mol
st = zeros(n,1); st(ii2) = -4*F; st(iNO2) = -1;
rate = react(j,C,params,op_cond);
if j ~= nj
    rateR = react(j+1,C,params,op_cond);
else
    rateR = rate;
end
w = 0.5;
gen  = st(i)*(w*rate+(1-w)*rateR)*dx/2;

if dt == 0
    acc = 0;
else
    if i == ii2
        a = params(6); Cdl = params(7);
        if j == nj
            dVdt = a*Cdl*((C(j,iv1)-C(j,iv2))-...
                (Cp(j,iv1)-Cp(j,iv2)))/dt;
        else
            dVdt = 0.5*a*Cdl*((C(j,iv1)-C(j,iv2))-...
                (Cp(j,iv1)-Cp(j,iv2)))/dt+...
                0.5*a*Cdl*((C(j+1,iv1)-C(j+1,iv2))-...
                (Cp(j+1,iv1)-Cp(j+1,iv2)))/dt;
        end
        acc = dVdt*dx/2;
    elseif i == iNO2
        R = 83.14; % cm3 bar / mol K
        T = op_cond(2); p = op_cond(5); CT = p/(T*R);
        if j == nj
            dcdt = CT*(C(j,iyO2)-Cp(j,iyO2))/dt;
        else
            dcdt = 0.5*CT*((C(j,iyO2)-Cp(j,iyO2))/dt)+...
                0.5*CT*((C(j+1,iyO2)-Cp(j+1,iyO2))/dt);
        end
        acc = dcdt*dx/2;
    end
end

NR = flux - gen + acc;
end