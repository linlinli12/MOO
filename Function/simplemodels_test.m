% Code name: parameter_screening_grid_2016.m
% Author: Thembi Mdlul
% Last update: Dec. 22, 2016
%
% This code load the sparse parameter grid for screening, which has a total
% of 2,00,000 LHS possible combninations of parameters.
% This function calls two other functions, the zebrafish model (FiniDiffMod2.m)
% and drosophila model (run_Drosophila_CORE).

% Model outputs for each parameter combination will be compared againest 
% experimental results at 4.7, 5.3,5.7,6.6 and 6.7 hpf for WT, CLF,
% 5.7 hpf NLF for zebrafish. 
% Model outputs for each parameter combination will be compared againest 
% experimental results for WT, Sof-/- and sog-/+ for drosophila . The total
% error sum of squares for each of these combinations will be calculated.
% Later we will set up a SSE threshold based on the initial values above,
% and figure out how many combinations (and which combinations) have SSE's
% <= the set threshold and have similar profiles to experiments.
function main
%% Drosophila Preliminaries 
clear;  clc;

addpath('./ExpData')
addpath('./Function')
%% Additional code
%%%%%%%%%% load Experimental results for zebrafish
load('ModeldataALL.mat');
% load('T57sym.mat');

%--------------------- Data at 4.7hpf-------------------------------------%
pWT47 = mean(wt47allout{1,3},2) - 9.1; % WT data minus background
pCLF47 = mean(chd47allout{1,3},2) -  9.1; % Chd mutant data minus background
% smooth data
% pWT47 = sgolayfilt(pWT47,10,35);
% pCLF47 = sgolayfilt(pCLF47,10,35);
% reflect data and average
pWT47 = (pWT47(1:18) + flipud(pWT47(19:36)))/2;
pCLF47 = (pCLF47(1:18) + flipud(pCLF47(19:36)))/2;
% variance for data
vWT47 = var(wt47allout{1,3} - 9.1,[],2);
temp_D = (vWT47(1:18) + flipud(vWT47(19:36)))/2;
vWT47 = sqrt(temp_D);
vCLF47 = var(chd47allout{1,3} - 9.1,[],2);
temp_D = (vCLF47(1:18) + flipud(vCLF47(19:36)))/2;
vCLF47 = sqrt(temp_D);

%---------------------- Data at 5.3hpf------------------------------------%
pWT53 = mean(wt53allout{1,3},2)-8.7; % WT data minus background
pCLF53 = mean(chd53allout{1,3},2)-8.7; % Chd mutant data minus background
% smooth data
% pWT53= sgolayfilt(pWT53,10,35);
% pCLF53 = sgolayfilt(pCLF53,10,35);
% reflect data and average
pWT53 = (pWT53(1:18) + flipud(pWT53(19:36)))/2;
pCLF53 = (pCLF53(1:18) + flipud(pCLF53(19:36)))/2;
% variance for data
vWT53 = var(wt53allout{1,3} - 8.7,[],2);
temp_D = (vWT53(1:18) + flipud(vWT53(19:36)))/2;
vWT53 = sqrt(temp_D);
vCLF53 = var(chd53allout{1,3} - 9.1,[],2);
temp_D = (vCLF53(1:18) + flipud(vCLF53(19:36)))/2;
vCLF53 = sqrt(temp_D);

%--------- Data at 5.7hpf-------------------------------------------------%
pWT57 = mean(wt57allout{1,3},2)-8.6;% WT data minus background
pCLF57 = mean(chd57allout{1,3},2)-8.6; % Chd mutant data minus background
pNLF57 = mean(NF57allout{1,3},2)-8.6;% Nog mutant data minus background
% smooth data
% pWT57= sgolayfilt(pWT57,10,35);
% pCLF57 = sgolayfilt(pCLF57,10,35);
% pNLF57 = sgolayfilt(pNLF57,10,35);

% reflect data and average
pWT57 = (pWT57(1:18) + flipud(pWT57(19:36)))/2;
pCLF57 = (pCLF57(1:18) + flipud(pCLF57(19:36)))/2;
pNLF57 = (pNLF57(1:18) + flipud(pNLF57(19:36)))/2;
% variance for data
vWT57 = var(wt57allout{1,3} - 8.6,[],2);
temp_D = (vWT57(1:18) + flipud(vWT57(19:36)))/2;
vWT57 = sqrt(temp_D);
vCLF57 = var(chd57allout{1,3} - 8.6,[],2);
temp_D = (vCLF57(1:18) + flipud(vCLF57(19:36)))/2;
vCLF57 = sqrt(temp_D);
vNLF57 = var(NF57allout{1,3} - 8.6,[],2);
temp_D = (vNLF57(1:18) + flipud(vNLF57(19:36)))/2;
vNLF57 = sqrt(temp_D);

%------------------ Data at 6.3hpf----------------------------------------%
pWT63 = mean(wt63allout{1,3},2)-9.5;% WT data minus background
pCLF63 = mean(chd63allout{1,3},2)-9.5; % Chd mutant data minus background
% smooth data
% pWT63= sgolayfilt(pWT63,10,35);
% pCLF63 = sgolayfilt(pCLF63,10,35);
% reflect data and average
pWT63 = (pWT63(1:18) + flipud(pWT63(19:36)))/2;
pCLF63 = (pCLF63(1:18) + flipud(pCLF63(19:36)))/2;
% variance for data
vWT63 = var(wt63allout{1,3} - 9.5,[],2);
temp_D = (vWT63(1:18) + flipud(vWT63(19:36)))/2;
vWT63 = sqrt(temp_D);
vCLF63 = var(chd63allout{1,3} - 9.5,[],2);
temp_D = (vCLF63(1:18) + flipud(vCLF63(19:36)))/2;
vCLF63 = sqrt(temp_D);

%--------------------- Data at 6.7hpf-------------------------------------%
pWT67 = mean(wt67allout{1,3},2)-10.1;% WT data minus background
pCLF67 = mean(chd67allout{1,3},2)-10.1; % Chd mutant data minus background
% smooth data
% pWT67= sgolayfilt(pWT67,10,35);
% pCLF67 = sgolayfilt(pCLF67,8,35);
% reflect data and average
pWT67 = (pWT67(1:18) + flipud(pWT67(19:36)))/2;
pCLF67 = (pCLF67(1:18) + flipud(pCLF67(19:36)))/2;
% variance for data
vWT67 = var(wt67allout{1,3} - 10.1,[],2);
temp_D = (vWT67(1:18) + flipud(vWT67(19:36)))/2;
vWT67 = sqrt(temp_D);
vCLF67 = var(chd67allout{1,3} - 10.1,[],2);
temp_D = (vCLF67(1:18) + flipud(vCLF67(19:36)))/2;
vCLF67 = sqrt(temp_D);


%Normalization constants

refmax2 = sort(pWT47,'descend');    % ref data set for scaling
refmax47 = mean(refmax2(1:5));

refmax2 = sort(pWT53,'descend');    % ref data set for scaling
refmax53 = mean(refmax2(1:5));

refmax2 = sort(pWT57,'descend');    % ref data set for scaling
refmax57 = mean(refmax2(1:5));

refmax2 = sort(pWT63,'descend');    % ref data set for scaling
refmax63 = mean(refmax2(1:5));

refmax2 = sort(pWT67,'descend');    % ref data set for scaling
refmax67 = mean(refmax2(1:5));


%%%%%%%%%% load Experimental results for Drosophila
load('chordinlike.mat'); % All chordilike sog data
temp{1} = load('WTData'); % All wildtype data
temp{2} = load('SogHomoData'); % All Sog homo data
temp{3} = load('SogHeterData'); % All Sog heter data

%% 
% data interpolation: use half embryo data
    xh = linspace(-1,0,18);
    dxh = linspace(-1,0,36);
    pWT47 = (interp1(xh,pWT47,dxh,'pchip'))/refmax47;
    pWT53 = (interp1(xh,pWT53,dxh,'pchip'))/refmax53;
    pWT57 = (interp1(xh,pWT57,dxh,'pchip'))/refmax57;
    pWT63 = (interp1(xh,pWT63,dxh,'pchip'))/refmax63;
    pWT67 = (interp1(xh,pWT67,dxh,'pchip'))/refmax67;
    pCLF47 = (interp1(xh,pCLF47,dxh,'pchip'))/refmax47;
    pCLF53 = (interp1(xh,pCLF53,dxh,'pchip'))/refmax53;
    pCLF57 = (interp1(xh,pCLF57,dxh,'pchip'))/refmax57;
    pCLF63 = (interp1(xh,pCLF63,dxh,'pchip'))/refmax63;
    pCLF67 = (interp1(xh,pCLF67,dxh,'pchip'))/refmax67;
    pNLF57 = (interp1(xh,pNLF57,dxh,'pchip'))/refmax57;
    % variance of data
    vWT47 = (interp1(xh,vWT47,dxh,'pchip'))/refmax47;
    vWT53 = (interp1(xh,vWT53,dxh,'pchip'))/refmax53;
    vWT57 = (interp1(xh,vWT57,dxh,'pchip'))/refmax57;
    vWT63 = (interp1(xh,vWT63,dxh,'pchip'))/refmax63;
    vWT67 = (interp1(xh,vWT67,dxh,'pchip'))/refmax67;
    vCLF47 = (interp1(xh,vCLF47,dxh,'pchip'))/refmax47;
    vCLF53 = (interp1(xh,vCLF53,dxh,'pchip'))/refmax53;
    vCLF57 = (interp1(xh,vCLF57,dxh,'pchip'))/refmax57;
    vCLF63 = (interp1(xh,vCLF63,dxh,'pchip'))/refmax63;
    vCLF67 = (interp1(xh,vCLF67,dxh,'pchip'))/refmax67;
    vNLF57 = (interp1(xh,vNLF57,dxh,'pchip'))/refmax57;

%% Zebrafish Preliminaries 
n = 36;   % number of nodes to evaluate finite difference equations
tRange = [0 4300 6500 7900 10100 11500];  % time interval to evaluate differential equations
  %%% geometry
zparameters.Ltot = 700;    % length of the embryo (1D assumption) microns
zparameters.Lven = 350;    % length of ventral region for Bmp  microns
zparameters.LvenXlr = 400; % length of ventral region for Tld  microns
zparameters.LdorC = 145;   % length of dorsal Chd expression region microns
zparameters.LdorN = 78;   % length of dorsal Nog expression region microns
  %%% Positive Feedback
zparameters.K = 0/60;         % Vmax for positive feedback (multiply by 0 to turn off)
zparameters.kM = 4000/60;       % kM parameter for Hill function
zparameters.nu = 2;       % cooperative parameter
  %%% Noggin fixed parameters  
load('./ExpData/parameters_zebra_noggin'); % Loard Noggin parameters previously optimized
zparameters.j3 = nondP(ind(1),30);  % production rate of Noggin       nM*s^-1
zparameters.k3 = nondP(ind(1),8);% binding rates for BMP ligand and Noggin          nM^-1*s^-1
zparameters.k_3 = nondP(ind(1),9); % unbinding rates for BMP ligand and Noggin          nM^-1*s^-1
zparameters.decN = nondP(ind(1),12);%0.15/60  * 10^(-2+4*p(7));    % decay rate of Nog             nM*m^-1
zparameters.decBN = nondP(ind(1),14);
zparameters.DN = nondP(ind(1),17);%10;       % diffusion rate of Noggin          (microns^2*s^-1)*60s*m^-1
zparameters.DBN = nondP(ind(1),19);%6;       % diffusion rate of [BMP,Nog]       (microns^2*s^-1)*60s*m^-1  
%%% Fixed Chd unbinding rate  
zparameters.k2 = nondP(ind(1),6); %binding rates for BMP ligand and Sog        nM^-1s^-1   (based on Mizutani 2005)
zparameters.tld_conc = 1;             % tld concentration           nM
%%% kinetic parameters 
zparameters.kendo1 = 0*1e-3; % endocytosis rate of [B,C]   m^-1  ##### kendo1 = kendo2
zparameters.kendo2 = 0*1e-3; % endocytosis rate of [BN,C]  m^-1  
zvar.parameters = zparameters;
zvar.tRange = tRange;
zvar.n = n; %number of nodes to evaluate finite difference equations
zvar.m = 2; %number of outputs
zvar.model = @zebrafish_model;     % model to run zebrafish
%% Drosophila Preliminaries 
dparameters.Ltot = 275;    % length of the embryo (1D assumption)       microns
dparameters.k2 = 0.095; %binding rates for BMP ligand and Sog        nM^-1s^-1   (based on Mizutani 2005)
dparameters.tld_conc = 1;             % tld concentration           nM
tRange = [0 3600];  % time interval to evaluate differential equations 5.7hpf
dvar.parameters = dparameters;
dvar.tRange = tRange;
dvar.n = n; %number of nodes to evaluate finite difference equations
dvar.m = 2; %number of outputs
dvar.model = @run_Drosophila_CORE;

%% Sampling uncertainty parameter space using LHS
% Define ranges
tic
d = 11;
znames = {'k_2','decB','decS','decBS','DB','DS','DBS','\lambda S',...
           '\lambda BS','j1','j2'};
prange = [0.0001 1; 0.00001 0.1; 0.00001 0.1; 0.00001 0.1;
    0.01 100; 0.01 100; 0.01 100; 0.0001 1; 0.0001 1
    0.01 100; 0.01 100];
prange = log10(prange);
sam_p = 1000000;  % number of samples
ALLP0 = lhsdesign(sam_p,d, 'criterion','maximin', 'criterion','correlation'); 
for k = 1:d
    ALLP0(:,k) = (prange(k,2)-prange(k,1))*ALLP0(:,k) + prange(k,1); 
end
ALLP0 = 10.^(ALLP0);
save('Parameters_screen','ALLP0');   % screening parameters are saved
toc
load('Parameters_screen');  
load('parameters_david')
kj = 1;
cnt_par = 1; 

%% Simulations of sampled parameters in badges of 2000
for i=1:250
    tic
    disp(['i = ', int2str(i)]);
   
%     %-------------- screening operations ----------------
    nreg = 2;
    save_p0 = ALLP0(cnt_par:cnt_par+nreg-1,:);    % select 2000 hundred parameters for each loop
%     save_p0 = parameters_new(cnt_par:cnt_par+nreg-1,:);  
    p0 = save_p0;
    save_p0 = p0;
    
%         p0 = [nondP(ind(1),7) nondP(ind(1),10) nondP(ind(1),11) nondP(ind(1),13)...
%         nondP(ind(1),15) nondP(ind(1),16) nondP(ind(1),18)...
%         nondP(ind(1),22) nondP(ind(1),23) nondP(ind(1),28) nondP(ind(1),29);...
%         0.1160 1.6408e-5 0.0037 5.0065e-4 0.1809 21.3097 7.4302 0.0021 0.3043 0.0113 7.3589];
%     save_p0 = p0;
    %%%%%%%%%% ODE solver
    %% ------------- Zebrafish Simulations
    %--- solve for WT------------------------------------------------------
    clear zc0;
    n = zvar.n; % number of time points
    m = zvar.m;  % number of outputs 
    zvar.tRange = [0 4300 6500 7900 10100 11500]; % simulation times
    ntp = numel(zvar.tRange)-1; % no. of simulation times
    zc0 = cell(1,2*n*m*ntp);
    
    [zc0{:}] = getStates(p0,zvar);
    
    zc0 = [zc0{:}];
    s = size(zc0,1);
    zc0 = reshape(zc0,s,ntp,[]);
    % time 4.7hpf
    ST1 = squeeze(zc0(:,1,:));
    B_WT1 = ST1(:,1:72); % BMP gradient
    C_WT1 = ST1(:,73:end);  % Chd gradient
    % time 5.3 hpf
    ST2 = squeeze(zc0(:,2,:));
    B_WT2 = ST2(:,1:72);
    C_WT2 = ST2(:,73:end);
    % time 5.7 hpf
    ST3 = squeeze(zc0(:,3,:));
    B_WT3 = ST3(:,1:72);
    C_WT3 = ST3(:,73:end); 
    % time 6.3 hpf
    ST4 = squeeze(zc0(:,4,:));
    B_WT4 = ST4(:,1:72);
    C_WT4 = ST4(:,73:end);
    % time 6.7 hpf
    ST5 = squeeze(zc0(:,5,:));
    B_WT5 = ST5(:,1:72);
    C_WT5 = ST5(:,73:end);
    
    %---- solve for CLF ---------------------------------------------------
      %%% production rates 
      clear zc0;
    p0 = save_p0;
    p0(:,11) = 0;  % production rate of Chordin      nM*m^-1
    ntp = numel(zvar.tRange)-1;
    zc0 = cell(1,2*n*m*ntp);
    [zc0{:}] = getStates(p0,zvar);
    zc0 = [zc0{:}];
    s = size(zc0,1);
    zc0 = reshape(zc0,s,ntp,[]);
    % time 4.7hpf
    ST1 = squeeze(zc0(:,1,:));
    B_CLF1 = ST1(:,1:72); % BMP gradient
    C_CLF1 = ST1(:,73:end); % Chd gradient
    % time 5.3 hpf
    ST2 = squeeze(zc0(:,2,:));
    B_CLF2 = ST2(:,1:72);
    C_CLF2 = ST2(:,73:end);
    % time 5.7 hpf
    ST3 = squeeze(zc0(:,3,:));
    B_CLF3 = ST3(:,1:72);
    C_CLF3 = ST3(:,73:end); 
    % time 6.3 hpf
    ST4 = squeeze(zc0(:,4,:));
    B_CLF4 = ST4(:,1:72);
    C_CLF4 = ST4(:,73:end);
    % time 6.7 hpf
    ST5 = squeeze(zc0(:,5,:));
    B_CLF5 = ST5(:,1:72);
    C_CLF5 = ST5(:,73:end);
    
     %---- solve for CLF heter ---------------------------------------------------
      %%% production rates 
      clear zc0;
    p0 = save_p0;
    p0(:,11) = 0.5*p0(:,11);  % production rate of Chordin      nM*m^-1
    zvar.tRange = [0 7900]; % simulation times
    ntp = numel(zvar.tRange)-1;
    zc0 = cell(1,2*n*m*ntp);
    [zc0{:}] = getStates(p0,zvar);
    zc0 = [zc0{:}];
    s = size(zc0,1);
    zc0 = reshape(zc0,s,ntp,[]);
    % time 5.7hpf
    ST1 = squeeze(zc0(:,1,:));
    B_CLFh1 = ST1(:,1:72); % BMP gradient
    C_CLFh1 = ST1(:,73:end); % Chd gradient

    %--- solve for NLF----------------------------------------------------
      %%% production rates 
      clear zc0;
    zparameters = zvar.parameters;
    zparameters.j3 = 0; % production rate of Noggin      nM*m^-1
    zvar.parameters = zparameters;
    zvar.tRange = [0 7900];
    p0 = save_p0;
    ntp = numel(zvar.tRange)-1;
    zc0 = cell(1,2*n*m*ntp);
    [zc0{:}] = getStates(p0,zvar);
    zc0 = [zc0{:}];
    B_NLF = zc0(:,1:72); % BMP gradient
    C_NLF = zc0(:,73:end); % Chd gradient

   %% ------------- Drosophila Simulations--------------------------------
    %------ Wild type simulations -----------------------------------------
      clear dc0;
    n = dvar.n; % number of time points
    m = dvar.m;  % number of outputs 
    ntp = numel(dvar.tRange)-1;
    dc0 = cell(1,2*n*m*ntp);
    p0 = save_p0;
    [dc0{:}] = getStates(p0,dvar);
    dc0 = [dc0{:}];
    dB_WT = dc0(:,1:72); % free BMP gradient
    dS_WT = dc0(:,73:144); % Sog gradient
    
    %---- solve for SLF ---------------------------------------------------
      clear dc0;
    p0 = save_p0;
    p0(:,11) = 0;  % production rate of Sog      nM*m^-1
    dc0 = cell(1,2*n*m*ntp);
    [dc0{:}] = getStates(p0,dvar);
    dc0 = [dc0{:}];
    dB_SHomo = dc0(:,1:72);
    dS_SHomo = dc0(:,73:144);
    
    %----- solve for SLF heter --------------------------------------------
      clear dc0;
    p0 = save_p0;
    p0(:,11) = 0.5* p0(:,11);  % production rate of Sog      nM*m^-1
    dc0 = cell(1,2*n*m*ntp);
    [dc0{:}] = getStates(p0,dvar);
    dc0 = [dc0{:}];
    dB_SHet = dc0(:,1:72);
    dS_SHet = dc0(:,73:144);

        %% Data fitting for zebrafish
    %%%----Save simulations--------
    % wild type
    Stored_Soln.zwt47 = B_WT1;  
    Stored_Soln.zwt53 = B_WT2;
    Stored_Soln.zwt57 = B_WT3;
    Stored_Soln.zwt63 = B_WT4;
    Stored_Soln.zwt67 = B_WT5;
    %
    Stored_Soln.zwtc47 = C_WT1;  
    Stored_Soln.zwtc53 = C_WT2;
    Stored_Soln.zwtc57 = C_WT3;
    Stored_Soln.zwtc63 = C_WT4;
    Stored_Soln.zwtc67 = C_WT5;
    % CLF 
    Stored_Soln.zclf47 = B_CLF1;  
    Stored_Soln.zclf53 = B_CLF2;
    Stored_Soln.zclf57 = B_CLF3;
    Stored_Soln.zclf63 = B_CLF4;
    Stored_Soln.zclf67 = B_CLF5;
    %
    Stored_Soln.zclfc47 = C_CLF1;  
    Stored_Soln.zclfc53 = C_CLF2;
    Stored_Soln.zclfc57 = C_CLF3;
    Stored_Soln.zclfc63 = C_CLF4;
    Stored_Soln.zclfc67 = C_CLF5;
    % CLF heter
    Stored_Soln.zhclf57 = B_CLFh1;
    Stored_Soln.zhclfc57 = C_CLFh1;
    % NLF
    Stored_Soln.znlf57 = B_NLF;
    Stored_Soln.znlfc57 = C_NLF;
    
         %%%%%%%%%% model results scaling
    B_CLF1   = B_CLF1(:,1:36).* repmat((1./(max(B_WT1(:,1:36),[],2))),1,36);
    B_CLF2   = B_CLF2(:,1:36).* repmat((1./(max(B_WT2(:,1:36),[],2))),1,36);
    B_CLF3   = B_CLF3(:,1:36).* repmat((1./(max(B_WT3(:,1:36),[],2))),1,36);
    B_CLF4   = B_CLF4(:,1:36).* repmat((1./(max(B_WT4(:,1:36),[],2))),1,36);
    B_CLF5   = B_CLF5(:,1:36).* repmat((1./(max(B_WT5(:,1:36),[],2))),1,36);
    B_hCLF  = B_CLFh1(:,1:36) .* repmat((1./(max(B_WT3(:,1:36),[],2))),1,36);
    B_NLF  = B_NLF(:,1:36) .* repmat((1./(max(B_WT3(:,1:36),[],2))),1,36);
    B_WT1   = B_WT1(:,1:36) .* repmat((1./(max(B_WT1(:,1:36),[],2))),1,36);
    B_WT2   = B_WT2(:,1:36) .* repmat((1./(max(B_WT2(:,1:36),[],2))),1,36);
    B_WT3   = B_WT3(:,1:36) .* repmat((1./(max(B_WT3(:,1:36),[],2))),1,36);
    B_WT4   = B_WT4(:,1:36) .* repmat((1./(max(B_WT4(:,1:36),[],2))),1,36);
    B_WT5   = B_WT5(:,1:36) .* repmat((1./(max(B_WT5(:,1:36),[],2))),1,36);

       
    %%%%%%%%%% Calculating sum of squares for zebrafish %%%%%%%%%%%%%%%%%
    SSE = zeros(nreg,17);
    % wild type
    SSE(:,1) = sum(((B_WT1 - repmat(pWT47,nreg,1))).^2,2);  
    SSE(:,2) = sum(((B_WT2 - repmat(pWT53,nreg,1))).^2,2); 
    SSE(:,3) = sum(((B_WT3 - repmat(pWT57,nreg,1))).^2,2);
    SSE(:,4) = sum(((B_WT4 - repmat(pWT63,nreg,1))).^2,2); 
    SSE(:,5) = sum(((B_WT5 - repmat(pWT67,nreg,1))).^2,2);  
    % CLF
    SSE(:,6) = sum(((B_CLF1 - repmat(pCLF47,nreg,1))).^2,2);
    SSE(:,7) = sum(((B_CLF2 - repmat(pCLF53,nreg,1))).^2,2);
    SSE(:,8) = sum(((B_CLF3 - repmat(pCLF57,nreg,1))).^2,2);
    SSE(:,9) = sum(((B_CLF4 - repmat(pCLF63,nreg,1))).^2,2);
    SSE(:,10) = sum(((B_CLF5 - repmat(pCLF67,nreg,1))).^2,2);
    % hCLF     
    SSE(:,11) = sum(((B_hCLF - repmat(pWT57,nreg,1))).^2,2);
    % NLF     
    SSE(:,12) = sum(((B_NLF - repmat(pNLF57,nreg,1))).^2,2);
    
    % the vector constant for linear approximation
    bd1 = abs(rdivide((B_WT1 - repmat(pWT47,nreg,1)),...
        repmat((diag(sqrt(36*((B_WT1 - repmat(pWT47,nreg,1))*...
        (B_WT1 - repmat(pWT47,nreg,1))')))),1,36)));
    bd2 = abs(rdivide((B_WT2 - repmat(pWT53,nreg,1)),...
        repmat((diag(sqrt(36*((B_WT2 - repmat(pWT53,nreg,1))*...
        (B_WT2 - repmat(pWT53,nreg,1))')))),1,36)));
    bd3 = abs(rdivide((B_WT3 - repmat(pWT57,nreg,1)),...
        repmat((diag(sqrt(36*((B_WT3 - repmat(pWT57,nreg,1))*...
        (B_WT3 - repmat(pWT57,nreg,1))')))),1,36)));
    bd4 = abs(rdivide((B_WT4 - repmat(pWT63,nreg,1)),...
        repmat((diag(sqrt(36*((B_WT4 - repmat(pWT63,nreg,1))*...
        (B_WT4 - repmat(pWT63,nreg,1))')))),1,36)));
    bd5 = abs(rdivide((B_WT5 - repmat(pWT67,nreg,1)),...
        repmat((diag(sqrt(36*((B_WT5 - repmat(pWT67,nreg,1))*...
        (B_WT5 - repmat(pWT67,nreg,1))')))),1,36)));
    
    bd6 = abs(rdivide((B_CLF1 - repmat(pCLF47,nreg,1)),...
        repmat((diag(sqrt(36*((B_CLF1 - repmat(pCLF47,nreg,1))*...
        (B_CLF1 - repmat(pCLF47,nreg,1))')))),1,36)));
    bd7 = abs(rdivide((B_CLF2 - repmat(pCLF53,nreg,1)),...
        repmat((diag(sqrt(36*((B_CLF2 - repmat(pCLF53,nreg,1))*...
        (B_CLF2 - repmat(pCLF53,nreg,1))')))),1,36)));
    bd8 = abs(rdivide((B_CLF3 - repmat(pCLF57,nreg,1)),...
        repmat((diag(sqrt(36*((B_CLF3 - repmat(pCLF57,nreg,1))*...
        (B_CLF3 - repmat(pCLF57,nreg,1))')))),1,36)));
    bd9 = abs(rdivide((B_CLF4 - repmat(pCLF63,nreg,1)),...
        repmat((diag(sqrt(36*((B_CLF4 - repmat(pCLF63,nreg,1))*...
        (B_CLF4 - repmat(pCLF63,nreg,1))')))),1,36)));
    bd10 = abs(rdivide((B_CLF5 - repmat(pCLF67,nreg,1)),...
        repmat((diag(sqrt(36*((B_CLF5 - repmat(pCLF67,nreg,1))*...
        (B_CLF5 - repmat(pCLF67,nreg,1))')))),1,36)));
    bd12 = abs(rdivide((B_NLF - repmat(pNLF57,nreg,1)),...
        repmat((diag(sqrt(36*((B_NLF - repmat(pNLF57,nreg,1))*...
        (B_NLF - repmat(pNLF57,nreg,1))')))),1,36)));

    
    Bd(:,1) = diag(bd1* diag(vWT47)*bd1');
    Bd(:,2) = diag(bd2* diag(vWT53)*bd2');
    Bd(:,3) = diag(bd3* diag(vWT57)*bd3');
    Bd(:,4) = diag(bd4* diag(vWT63)*bd4');
    Bd(:,5) = diag(bd5* diag(vWT67)*bd5');
    
    Bd(:,6) = diag(bd6* diag(vCLF47)*bd6');
    Bd(:,7) = diag(bd7* diag(vCLF53)*bd7');
    Bd(:,8) = diag(bd8* diag(vCLF57)*bd8');
    Bd(:,9) = diag(bd9* diag(vCLF63)*bd9');
    Bd(:,10) = diag(bd10* diag(vCLF67)*bd10');
    
    Bd(:,11) = diag(bd3* diag(vWT57)*bd3');
    Bd(:,12) = diag(bd12* diag(vNLF57)*bd12');
   %% Data fitting for drosophila
   %%%------Save simulations----- 
    Stored_Soln.dwt = dB_WT;
    Stored_Soln.dwts = dS_WT;
    Stored_Soln.dshomo = dB_SHomo;
    Stored_Soln.dshet = dB_SHet;
    
    dB_WT = dB_WT(:,1:36);     % BMP signal drosophila
    dB_SHomo = dB_SHomo(:,1:36);     % BMP signal drosophila
    dB_SHet = dB_SHet(:,1:36);     % BMP signal drosophila
    
    % temp is a temporary structure for drosophila data
    % wild type 
    dxh = linspace(min(temp{1}.data.x),max(temp{1}.data.x),100);
    dWT = interp1(temp{1}.data.x,temp{1}.data.y,dxh,'pchip');
    xh = linspace(-1,0,36);
    dWT = interp1(dxh,dWT,xh,'pchip');
    dWT = repmat(dWT,nreg,1);
    dB_WTt = dB_WT./repmat(max(dB_WT,[],2),1,36);
    % variance 
    dst1 = interp1(temp{1}.data.x,temp{1}.data.st,xh,'pchip');
    % apply weights for sog mutant
    wgts = ones(1,numel(xh));
    wgts(xh > -0.6 & xh < -0.45) = 0.1;
    
    % sog loss of function
    dxh = linspace(min(temp{2}.data.x),max(temp{2}.data.x),100);
    dSLF = interp1(temp{2}.data.x,temp{2}.data.y,dxh,'pchip');
    dSLF = interp1(dxh,dSLF,xh,'pchip');
    dst2 = 0.1*dSLF;
    dSLF = repmat(dSLF,nreg,1);
    dB_SHomot = dB_SHomo./repmat(max(dB_WT,[],2),1,36);
    % variance 
    dst2(xh > -0.6 & xh < -0.45) = 0.1;
    dst2(dst2 == 0) = dst1(dst2 == 0);
    
    % sog loss of function heter
    dxh = linspace(min(temp{3}.data.x),max(temp{3}.data.x),100);
    dSHet = interp1(temp{3}.data.x,temp{3}.data.y,dxh,'pchip');
    dSHet = interp1(dxh,dSHet,xh,'pchip');
    dSHet = repmat(dSHet,nreg,1);
    dB_SHett = dB_SHet./repmat(max(dB_WT,[],2),1,36);
    % variance 
    dst3 = interp1(temp{3}.data.x,temp{3}.data.st,xh,'pchip');


      %-------------with chordinlike data--------------------------
    
    % sog loss of function
    dxh = linspace(min(temp{2}.data.x),max(temp{2}.data.x),100);
    dSLFt = interp1(temp{2}.data.x,temp{2}.data.y,dxh,'pchip');
    dSLFt = interp1(dxh,dSLFt,xh,'pchip');
    SHomot = dB_SHomo./repmat(max(dB_WT,[],2),1,36);
    dSLFt = dSLFt/max(chordinlike.y);
    dSLFt = repmat(dSLFt,nreg,1);
    % variance 
    dst5 = dst2*(max(dSLFt(1,:),[],2)/max(dSLF(1,:),[],2));
%     dst5 = interp1(temp{2}.data.x,temp{2}.data.st,xh,'pchip');
%     ind = find(dst5);       dst4(dst5 == 0) = 0.1*dst5(ind(1)); 
    
    % Peluso WT
    dxh = linspace(min(chordinlike.x),max(chordinlike.x),1000);
    dWTt = interp1(chordinlike.x,chordinlike.y,dxh,'pchip');
    dWTt = dWTt/max(dWTt);    
    xh = linspace(-1,0,36);
    dWTt = interp1(dxh,dWTt,xh,'pchip'); 
    dWTt = repmat(dWTt,nreg,1);
    WTt = dB_WT./repmat(max(dB_WT,[],2),1,36);
    % variance 
    dst4 = interp1(chordinlike.x,chordinlike.st,xh,'pchip');
    ind = find(dst4);       dst4(dst4 == 0) = 0.1*dst4(ind(1)); 
    
    
    %%%%%%%%%% Calculating sum of squares for zebrafish %%%%%%%%%%%%%%%%%
    SSE(:,13) = sum((((dB_WTt - dWT)).^2),2); 
    SSE(:,14) = (((dB_SHomot - dSLF)).^2)*wgts';
    SSE(:,15) = sum((((dB_SHett - dSHet)).^2),2);

    % chordinlike fitness
    SSE(:,16) = sum((((WTt - dWTt)).^2),2);  
    SSE(:,17) = (((SHomot - dSLFt)).^2)*wgts'; 
    SSE = sqrt((1/36)*SSE);
    
    % the vector constant for linear approximation
    bd13 = abs(rdivide((dB_WTt - dWT),repmat((diag(sqrt(36*((dB_WTt - dWT)*...
        (dB_WTt - dWT)')))),1,36)));
    bd14 = abs(rdivide((dB_SHomot - dSLF),repmat((diag(sqrt(36*((dB_SHomot - dSLF)*...
        (dB_SHomot - dSLF)')))),1,36)));
    bd15 = rdivide((dB_SHett - dSHet),repmat((diag(sqrt(36*(dB_SHett - dSHet)*...
        (dB_SHett - dSHet)'))),1,36));
    bd16 = abs(rdivide((WTt - dWTt),repmat((diag(sqrt(36*((WTt - dWTt)*...
        (WTt - dWTt)')))),1,36)));
    bd17 = abs(rdivide((SHomot - dSLFt),repmat((diag(sqrt(36*((SHomot - dSLFt)*...
        (SHomot - dSLFt)')))),1,36)));
    
    Bd(:,13) = diag(bd13* diag(dst1.^2) *bd13');
    Bd(:,14) = diag(bd14* diag(dst2.^2) *bd14');
    Bd(:,15) = diag(bd15* diag(dst3.^2) *bd15');
    Bd(:,16) = diag(bd16* diag(dst4.^2) *bd16');
    Bd(:,17) = diag(bd17* diag(dst5.^2) *bd17');
    
    % save
    parameters_store = save_p0;
    fname=sprintf('SSE_scenario_EXTENDEDMODELS_%d.mat',kj);
    save(fname,'SSE','parameters_store','Stored_Soln','Bd');
    kj=kj+1;
    cnt_par = cnt_par+nreg;
    toc
end


end
