%% Main Algorithm, computing the option values and optimal exercise boundaries using FDM+PDAS

function [V1,V2,b_1,b_2] = FDM_PDAS(M,N,dx,dt,x,t,L,sigma,r,a,x0,t0,K)
    format long

    r1 = r(1);
    r2 = r(2);
    sigma1 = sigma(1);
    sigma2 = sigma(2);
    beta1 = 1/2-r1/(sigma1^2);
    beta2 = 1/2-r2/(sigma2^2);
    alpha1 = -(r1-a(1,1))-1/(2*sigma1^2)*(r1-sigma1^2/2)^2;
    alpha2 = -(r2-a(2,2))-1/(2*sigma2^2)*(r2-sigma2^2/2)^2;
    xi = alpha2-alpha1;
    eta = beta2-beta1;
    alpha = dt/(dx^2);  %the net ratio=alpha
    
    % construct matrix D  

    B1_main_diag =  1+alpha*sigma1^2*ones(N-1,1);
    B2_main_diag =  1+alpha*sigma2^2*ones(N-1,1);
    B1_sub_diag = -alpha*sigma1^2/2*ones(N-2,1);
    B2_sub_diag = -alpha*sigma2^2/2*ones(N-2,1);
    B1 = diag(B1_main_diag, 0) + diag(B1_sub_diag, -1) + diag(B1_sub_diag, 1);
    B2 = diag(B2_main_diag, 0) + diag(B2_sub_diag, -1) + diag(B2_sub_diag, 1);
    
    D = [B1 zeros(size(B1));zeros(size(B1)) B2];  
    
    % construct matrice A1 and A2
    A1 = -dt*a(1,2)*exp(xi*t0+eta*x0)*(exp((1:N-1)'*dx*eta));
    A2 = -dt*a(2,1)*exp(-(xi*t0+eta*x0))*(exp(-(1:N-1)'*dx*eta));

    % Initialize W and b for option values and optimal exercise boundaries
    W1 = zeros(M+1,N+1);
    W2 = zeros(M+1,N+1);
    b_1 = ones(M+1, 1)*K;
    b_2 = ones(M+1, 1)*K;
    
    S = exp(x);

    
    % Assign initial and boundary values to W 
    W1(1,:) = exp(-beta1*x).*max(K-S,0);
    W2(1,:) = exp(-beta2*x).*max(K-S,0);
    W1(:,N+1) = 0;
    W2(:,N+1) = 0;
    W1(:,1) = exp(-(alpha1*t+beta1*(-L))).*max(K-exp(-L),0);
    W2(:,1) = exp(-(alpha2*t+beta2*(-L))).*max(K-exp(-L),0);

    F0_1 = -alpha*sigma1^2/2*exp(-(alpha1*t0+beta1*(-L)))*max(K-exp(-L),0);
    F0_2 = -alpha*sigma2^2/2*exp(-(alpha2*t0+beta2*(-L)))*max(K-exp(-L),0);
    G0_1 = (exp(-beta1*x(2:N)).*max(K- S(2:N),0))';
    G0_2 = (exp(-beta2*x(2:N)).*max(K- S(2:N),0))';
    
    % time iterations
    for i = 2:M+1
        U1 = W1(i-1,2:N);
        U2 = W2(i-1,2:N);
        F = zeros(2*N-2,1);
        F(1) = F0_1*exp(-alpha1*(i-1)*dt);  
        F(N) = F0_2*exp(-alpha2*(i-1)*dt);        
        F(1:N-1,1)=F(1:N-1,1)-U1'+exp(xi*dt*(i-2))*A1.*(U2');
        F(N:end,1)=F(N:end,1)+exp(-xi*dt*(i-2))*A2.*(U1')-U2';

        G1 = exp(-alpha1*(i-1)*dt)*G0_1;
        G2 = exp(-alpha2*(i-1)*dt)*G0_2;
        G = [G1;G2];
        lambda=max(D*G+F,zeros(size(F)));
        [Phi_end,b1,b2] = PDAS(D,F,G,N,lambda);  %PDAS iteration
        W1(i,2:N) = Phi_end(1:N-1);
        W2(i,2:N) = Phi_end(N:2*N-2);
        b_1(i)=S(b1);
        b_2(i)=S(b2);
    end
    % store the optimal exercise boundaries
    b_1 = b_1(M+1:-1:1);
    b_2 = b_2(M+1:-1:1);  
    
    % Transform from W to V, the option values; see Eq. (6)
    V1 = exp(beta1*ones(M+1,1)*x+alpha1*t'*ones(1,N+1)).*W1;      
    V2 = exp(beta2*ones(M+1,1)*x+alpha2*t'*ones(1,N+1)).*W2;  
end
