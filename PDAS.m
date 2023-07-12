%% PDAS algorithm
% See Algorithm 1

function [Phi_end,b1,b2] = PDAS(D,F,G,Nx,lambda)  

format long
    rho = 1e8;
    err = 1e-8;
    Phi_pre = zeros(size(F));   
    Phi_cur = G;   
    k = 1;
    while norm(Phi_pre-Phi_cur,inf)>err
        Phi_pre = Phi_cur;
        IS = find((lambda+rho*(G-Phi_pre))<=0);
        AS = find((lambda+rho*(G-Phi_pre))>0);
        Phi_cur(AS) = G(AS);
        lambda(IS) = zeros(size(IS));
        Phi_cur(IS) = D(IS,IS)\(lambda(IS)-F(IS)-D(IS,AS)*G(AS));
        lambda(AS) = D(AS,:)*Phi_cur+F(AS);
        k = k+1; 
    end  
    b1 = min(IS);
    b2 = max(AS)-(Nx-1)+1;
    Phi_end = Phi_cur;
end

