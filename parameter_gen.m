close all
clear all
% this code is for crate parameter set with in the individual parameter
% space

folder_n = 200; % folders for paralle compution
case_n = 5000;% simulaiton cased in each folder
file_name = 'ParaSet_base_200x5000.mat';

rng('shuffle')
tic
d = 11;
znames = {'k2','decB','decS','decBS','DB','DS','DBS','lambdaS',...
           'lambdaBS','j1','j2'};
% setup individual parameter range
prange = [0.0001 1; 0.00001 0.1; 0.00001 0.1; 0.00001 0.1;
    0.01 100; 0.01 100; 0.01 100; 0.0001 1; 0.0001 1
    0.01 100; 0.01 100];
prange = log10(prange);
sam_p = folder_n*case_n;  % number of samples

% best parameter for simple model
% latin hyper cube sampling 
screen_param = lhsdesign(sam_p,d, 'criterion','maximin', 'criterion','correlation'); 
% 
Z_ALLP0 = zeros(sam_p,d);
for k = 1:11
    % sample parameter in log space
    p_log = logspace(prange(k,1),prange(k,2),sam_p);
    Z_ALLP0(:,k) =p_log(round(screen_param(:,k).*sam_p)); 
end
save(file_name,'Z_ALLP0')
