function [BR, S, BS, B, X, T] = run_Drosophila_CORE_feedback(n, tRange, parameters,pars)
% function name: FiniDiffMod2.m
% author: Wei Dou
%
% David Umulis Updates on 07/01/2014
% Last Update: 12/3/2013
%
% functionality:
%    This function uses a 1D finite difference method to solve the ODE for
%    Zebrafish BMP gradient formation during embryogenesis. 
%    Based on the primary function FiniteTest.m for the evaluation of differential 
%    equations using the 4th order runge-kutta technique with a fixed step-size.
%
% ODE system:
%
%    d[B]/dt = DB*(d^2[B]/dx^2)
%            - kon2*[B]*[C] + koff2*[BC] 
%            - kon3*[B]*[N] + koff3*[BN] 
%
%            - decB*[B] + etaB + lambdaBC*Tld*[BC]             
%
%    d[C]/dt = DC*(d^2[C]/dx^2)
%            - kon2*[B]*[C] + koff2*[BC] 
%            - decC*[C] + etaC - lambdaC*Tld*[C]
%
%    d[N]/dt = DN*(d^2[N]/dx^2) 
%            - kon3*[B]*[N] + koff3*[BN]
%            - decN*[N] + etaN
%
%    d[BC]/dt = DBC*(d^2[BC]/dx^2)
%             + kon2*[B]*[C] - koff2*[BC] 
%             - lambdaBC*Tld*[BC]
%
%    d[BN]/dt = DBN*(d^2[BN]/dx^2)
%             + kon3*[B]*[N] - koff3*[BN]
%
% function input:
%    n - number of nodes to evaluate finite difference equations
%    tRange - time interval to evaluate differential equations, eg. [0 5000]
%    parameters - a stucture storing input parameter values like diffusion
%                 rates, and reaction rates, etc. 
%        .Ltot - length of the embryo (1D assumption)       microns
%        .Lven - length of ventral region (1D assumption)   microns
%        .LdorC - length of dorsal Chd expression region    microns
%        .LdorN - length of dorsal Nog expression region    microns
%        .DB - diffusion rate of ligand (BMP)    (microns^2*s^-1)*60s*m^-1
%        .DC - diffusion rate of Chordin       (microns^2*s^-1)*60s*m^-1
%        .DN - diffusion rate of Noggin          (microns^2*s^-1)*60s*m^-1
%        .DBC - diffusion rate of [BMP,Chd]    (microns^2*s^-1)*60s*m^-1
%        .DBN - diffusion rate of [BMP,Nog]      (microns^2*s^-1)*60s*m^-1
%        .k2 - binding rates for BMP ligand and Chordin      nM^-1*m^-1
%        .k_2 - unbinding rates for BMP ligand and Chordin   m^-1
%        .k3 - binding rates for BMP ligand and Noggin       nM^-1*m^-1
%        .k_3 - unbinding rates for BMP ligand and Noggin    m^-1
%        .k4 - binding rates for BMP ligand and surface component  nM^-1*m^-1
%        .k_4 - unbinding rates for BMP ligand and surface componetn  m^-1 
%        .k6 - binding rates for [BMP,Nog] and surface component  nM^-1*m^-1
%        .k_6 - unbinding rates for [BMP,Nog] and surface component  m^-1
%        .kendo1 - endocytosis rate of [B,C]   m^-1
%        .kendo2 - endocytosis rate of [BN,C]  m^-1
%        .decB - decay rate of Ligand (BMP)    m^-1  
%        .decC - decay rate of Chd             m^-1  
%        .decN - decay rate of Nog             m^-1  
%        .j1 - production rate of BMP          nM*m^-1  
%        .j2 - production rate of Chordin      nM*m^-1
%        .j3 - production rate of Noggin       nM*m^-1
%        .lambda_tld_C - tld processing rate of Chd     nM^-1*m^-1
%        .lambda_tld_BC - tld processing rate of BChd   nM^-1*m^-1
%        .tld_conc - tld concentration         nM
%        .Stot - total surface component       nM

% NOTE: make n that can exactly divided by Ltot for easy calculation
%
% function output:
%    B - vector storing distribution of BMP 
%    S - vector storing distribution of Sog
%    BS - vector storing distribution of BMP-Sog complex
%    X - vector storing postional information (micron)
%    T - vector storing time information


fdAdt = @dAdt;      %Create a function handle for faster evaluation of the differential equation
%-------------Varying parameters
parameters.k2 = pars(1);
parameters.k_2 = pars(1);%0.28e-3 * 10^(-2+p(1));        % binding rates for BMP ligand and Chordin          nM^-1*s^-1 
parameters.decB = pars(2);%0.005/60 * 10^(-2+p(5));   % decay rate of Bmp             nM*m^-1
parameters.decS = pars(3);%0.01/60 * 10^(-2+p(6));       % decay rate of Chd             nM*m^-1 
parameters.decBS = pars(4);  
parameters.DB = pars(5);%30;       % diffusion rate of ligand (BMP)    (microns^2*s^-1)*60s*m^-1
parameters.DS = pars(6);%20;       % diffusion rate of Chordin         (microns^2*s^-1)*60s*m^-1
parameters.DBS = pars(7);%13;      % diffusion rate of [BMP,Chd]       (microns^2*s^-1)*60s*m^-1    
  %% Tolloid behavior
parameters.lambda_tld_S = pars(8);%     % tld processing rate of Chd  nM^-1*m^-1
parameters.lambda_tld_BS = pars(9);%5/60;   % tld processing rate of LC   nM^-1*m^-1

%%% production rates
parameters.j1 = pars(10);        % production rate of BMP          nM*m^-1 
parameters.j2 = pars(11);  % production rate of Chordin      nM*m^-1

%--------------Initial condition vectors---------------------------
B0 = zeros(1,n);    %Dpp
S0 = zeros(1,n);    %Initialize vectors for Sog
BS0 = zeros(1,n);    %Initialize vectors for Sog/Tsg/Dpp
BR0 = zeros(1,n);    % Initialize vectors for receptors
C0 = zeros(1,n);    % Initialize vectors for Cv2
BC0 = zeros(1,n);    % Initialize vectors for Cv2/Dpp
BCR0 = zeros(1,n);    % Initialize vectors for Cv2/Dpp/receptors

initial = [B0 S0 BS0 BR0 C0 BC0 BCR0];  % Initial condition vector

%----------Solve ODE------------------------
options = odeset('RelTol',1e-9);
[T, D] = ode15s(fdAdt,tRange,initial,options,n,parameters);    

%------------------------- store data in vectors --------------------------
B1 = D(:, 1:n);      
S1 = D(:, n+1:2*n);
BS1 = D(:, 2*n+1:3*n);
BR1 = D(:, 3*n+1:4*n);
C1 = D(:, 4*n+1:5*n);
BC1 = D(:, 5*n+1:6*n);
BCR1 = D(:, 6*n+1:7*n);

%------------------------- Reflect the solution ---------------------------
B2 = zeros(size(B1));
S2 = zeros(size(S1));
BS2 = zeros(size(BS1));
BR2 = zeros(size(BR1));
C2 = zeros(size(C1));
BC2 = zeros(size(BC1));
BCR2 = zeros(size(BCR1));

for i=1:n
    B2(:,i) = B1(:,n+1-i);
    S2(:,i) = S1(:,n+1-i);
    BS2(:,i) = BS1(:,n+1-i);
    BR2(:,i) = BR1(:,n+1-i);
    C2(:,i) = C1(:,n+1-i);
    BC2(:,i) = BC1(:,n+1-i);
    BCR2(:,i) = BCR1(:,n+1-i);
end

B = [B1 B2];
S = [S1 S2];
BS = [BS1 BS2];
BR = [BR1 BR2];
C = [C1 C2];
BC = [BC1 BC2];
BCR = [BCR1 BCR2];

%--------------------------- make x data vector --------------------------- 
start = -parameters.Ltot + parameters.Ltot/n;
X = start:(parameters.Ltot/n):parameters.Ltot;
 
end

%--------------------Differential Equation Vector -------------------------
function dY = dAdt(t,Y,n,param)
% This function is the set of differential equations.  
  %%% geometry
Ltot = param.Ltot;   % length of the embryo (1D assumption)       microns
L2 = 150;           %microns from VM for dorsal BMP boundary
L3 = 50;            %microns from VM for ventral Sog boundary

n2 = round(L2*n/Ltot);      %275 microns for circumference/2 estimate
n3 = round(L3*n/Ltot);

%%%
k2 = param.k2;                %Sog?Tsg bingding bmp    (nM*min)^-1
k02 = param.k_2;               %reverse rxn rate    (min)^-1
DB =  param.DB;             %Diffusion coefficient for B (microns^2/sec)
DS =  param.DS;               %Diffusion coefficient for S 
DBS = param.DBS;               %Diffusion coefficient for BS
phi_B = param.j1;
phi_Sog = param.j2;
lambda_BS =  param.lambda_tld_BS;           %tld processing rate for BS(nM*min)^-1
lambda_S =  param.lambda_tld_S;           %tld processing rate for S (nM*min)^-1
decB = param.decB;              %decay rate B
decS = param.decS;              %decay rate S
decBS = param.decBS;              %decay rate BS
tld_conc = param.tld_conc;             %tld conc. nM
cop_r = 2;
k4 = 0.0008;                %Sog?Tsg bingding bmp    (nM*min)^-1
k_4 = 0.1333;               %reverse rxn rate    (min)^-1
decE = 0.0005;
Rtot = 394.3;
lamda_r = 0.5225;
Kh = 61.83;
k5 = 0.03333;
k_5 = 0.066667;
k6 = 0.0166667;
k_6 = 0.3333;
k7 = 0.00416667;
k_7 = 0.33333;

dx = 275/n;               %Distance between nodes
etaSog = zeros(1,n);      %Initialize prepatterns for secretion and Tld
etaB = zeros(1,n);
Tld = zeros(1,n);
% etaTsg = zeros(1,n);

for i = 1:n
    etaSog(i) = phi_Sog*double((i>(n3-1))*(i<(n2+1)));
    etaB(i) = phi_B*double((i>(n2-1)));
    Tld(i) = tld_conc*double((i>(n2-1)));
end

%%
B = Y(1:n);
S = Y(n+1:2*n);
BS = Y(2*n+1:3*n);
BR = Y(3*n+1:4*n);
C = Y(4*n+1:5*n);
BC = Y(5*n+1:6*n);
BCR = Y(6*n+1:7*n);

dB = zeros(1,n);
dS = zeros(1,n);
dBS = zeros(1,n);
dBR = zeros(1,n);
dC = zeros(1,n);
dBC = zeros(1,n);
dBCR = zeros(1,n);



    R = Rtot - BR(1)- BCR(1);
    
    dB(1) = DB/dx^2*[B(1)-2*B(1)+B(2)] + etaB(1) - k2*B(1)*S(1)...
            + k02*BS(1) + lambda_BS*Tld(1)*BS(1) - decB*B(1)...
            - k4*B(1)*R + k_4*BR(1) - k5*B(1)*C(1) + k_5*BC(1);   %
        
    dS(1) = DS/dx^2*[S(1)-2*S(1)+S(2)] - k2*B(1)*S(1) + k02*BS(1)...
            + etaSog(1) - lambda_S*Tld(1)*S(1) - decS*S(1);
    
    dBS(1) = DBS/dx^2*[BS(1)-2*BS(1)+BS(2)] + k2*B(1)*S(1) - k02*BS(1)...
            - lambda_BS*Tld(1)*BS(1) - decBS*BS(1);
        
    dC(1) = ((lamda_r*(BR(1))^(cop_r))/(Kh^(cop_r) + (BR(1))^(cop_r)))...
        - k5*B(1)*C(1) + k_5*BC(1) - k6*BR(1)*C(1) + k_6*BCR(1) - decE*C(1);
    
    dBC(1) = k5*B(1)*C(1) - k_5*BC(1) - k6*BC(1)*R + k_6*BCR(1) - decE*BC(1);
    
    dBCR(1) = k6*BC(1)*R - k_6*BCR(1) + k7*BR(1)*C(1) - k_7*BCR(1) - decE*BCR(1);
    
    dBR(1) = k4*B(1)*R - k_4*BR(1) + k_7*BCR(1) - k7*BR(1)*C(1) - decE*BR(1);
    
    for i=2:1:n-1
        
        R = Rtot - BR(i)- BCR(i);
    
        dB(i) = DB/dx^2*[B(i-1)-2*B(i)+B(i+1)] + etaB(i) - k2*B(i)*S(i)...
                + k02*BS(i) + lambda_BS*Tld(i)*BS(i) - decB*B(i)...
                - k4*B(i)*R + k_4*BR(i) - k5*B(i)*C(i) + k_5*BC(i);   %

        dS(i) = DS/dx^2*[S(i-1)-2*S(i)+S(i+1)] - k2*B(i)*S(i) + k02*BS(i)...
                + etaSog(i) - lambda_S*Tld(i)*S(i) - decS*S(i);

        dBS(i) = DBS/dx^2*[BS(i-1)-2*BS(i)+BS(i+1)] + k2*B(i)*S(i) - k02*BS(i)...
                - lambda_BS*Tld(i)*BS(i) - decBS*BS(i);

        dC(i) = ((lamda_r*(BR(i))^(cop_r))/(Kh^(cop_r) + (BR(i))^(cop_r)))...
            - k5*B(i)*C(i) + k_5*BC(i) - k6*BR(i)*C(i) + k_6*BCR(i) - decE*C(i);

        dBC(i) = k5*B(i)*C(i) - k_5*BC(i) - k6*BC(i)*R + k_6*BCR(i) - decE*BC(i);

        dBCR(i) = k6*BC(i)*R - k_6*BCR(i) + k7*BR(i)*C(i) - k_7*BCR(i) - decE*BCR(i);
    
        dBR(i) = k4*B(i)*R - k_4*BR(i) + k_7*BCR(i) - k7*BR(i)*C(i) - decE*BR(i);
        
    end
    
    R = Rtot - BR(n)- BCR(n);
    
    dB(n) = DB/dx^2*[B(n-1)-2*B(n)+B(n)] + etaB(n) - k2*B(n)*S(n)...
            + k02*BS(n) + lambda_BS*Tld(n)*BS(n) - decB*B(n)...
            - k4*B(n)*R + k_4*BR(n) - k5*B(n)*C(n) + k_5*BC(n);   %
        
    dS(n) = DS/dx^2*[S(n-1)-2*S(n)+S(n)] - k2*B(n)*S(n) + k02*BS(n)...
            + etaSog(n) - lambda_S*Tld(n)*S(n) - decS*S(n);
    
    dBS(n) = DBS/dx^2*[BS(n-1)-2*BS(n)+BS(n)] + k2*B(n)*S(n) - k02*BS(n)...
            - lambda_BS*Tld(n)*BS(n) - decBS*BS(n);
        
    dC(n) = ((lamda_r*(BR(n))^(cop_r))/(Kh^(cop_r) + (BR(n))^(cop_r)))...
        - k5*B(n)*C(n) + k_5*BC(n) - k6*BR(n)*C(n) + k_6*BCR(n) - decE*C(n);
    
    dBC(n) = k5*B(n)*C(n) - k_5*BC(n) - k6*BC(n)*R + k_6*BCR(n) - decE*BC(n);
    
    dBCR(n) = k6*BC(n)*R - k_6*BCR(n) + k7*BR(n)*C(n) - k_7*BCR(n) - decE*BCR(n);
    
    dBR(n) = k4*B(n)*R - k_4*BR(n) + k_7*BCR(n) - k7*BR(n)*C(n) - decE*BR(n);

 dY = [dB dS dBS dBR dC dBC dBCR]';

end