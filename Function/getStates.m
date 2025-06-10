function varargout = getStates(Points,var)
%GETSTATES Computes model states at specified time points
%   This function computes model states at specified spatial points (tp) for
%   the parameters (p) of the specified model (params).  This function is
%   designed to use the Parallel Toolbo with vectorized inputs
%   and multiple outputs, and for parallel computing using 'matlabpool'.
%   
%   Inputs:
%   p = an M-by-d matrix of parameter values.
%   var = structure containing handle for model file (e.g. @modelODEs),
%         vector of initial conditions, etc.
%   
%   Outputs:
%   varargout = cell structure, each cell contains state trajectories at a
%       specified time point.

% Extract working variables
p = Points;

% Extract working variables
model = var.model;
n = var.n;
m = var.m;
parameters = var.parameters;
tRange = var.tRange;
ntp = (numel(tRange)-1);
numOutputs = 2*n*m*ntp;
x2 = zeros(size(p,1),numOutputs);

% Compute model states 
parfor i = 1:size(p,1)
    [B,C] = model(n, tRange, parameters,p(i,:));

    D = [B C];
    if ntp == 1
        xt = D(end,:); x = xt(:)';
    else
        xt = D(2:end,:); x = xt(:)';
    end
    x2(i,:) = x;
end
% Define output variable
g = cell(1,numOutputs);
for k = 1:numOutputs, g{k} = x2(:,k); end
varargout = g;

end