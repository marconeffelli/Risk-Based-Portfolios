function weights = RiskPortfolios(S,type)

% function RiskPortfolios
% 
%
% Description:
% this function calculates the optimal weights of four risk-based
% portfolios, namely: the Minimum Variance, the Inverse Volatility, the
% Equal-Risk-Contribution and the Maximum Diversification.
% Weights are subject to the long-only constraints
%
% Inputs:
% - a covariance matrix of asset returns, S
% - the portfolio type: MV, IV, ERC, MD
%
% Usage:
% weights = RiskPortfolios(S,'MV')
%
% This version: 10/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is released under the BSD 2-clause license.

% Copyright (c) 2018, Marco Neffelli 
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%to improve -- add check for S is squared
p = size(S,1);

% Minimum Variance (MV)
if type == MV
m=zeros(p,1);
port = Portfolio('assetmean', m, 'assetcovar', S, ...
'lowerbudget', 1, 'upperbudget', 1, 'lowerbound', 0);
weights = estimateFrontier(port,1);
end

% Equal Risk Budget or Inverse volatility weighted portfolio (IV)
if type == IV
std_all = diag(sqrt(S));
inv_std_all = 1./std_all;
sum_inv_std_all = sum(inv_std_all);
weights = inv_std_all/sum_inv_std_all;
end

% Equally Risk Contribution (ERC)
if type == ERC
LB = zeros(1,p);
UB = ones(1,p);
beq =1;
f1 = @(w1) ERC_ObjectiveFunction(w1,S);  
w0 = 1/p*ones(p,1);
weights = fmincon(f1,w0,[],[],UB,beq,LB,UB',[]);
end

% Most Diversified (MD)
if type == MDP
LB = zeros(1,p);
UB = ones(1,p);
beq =1;
std_all = diag(sqrt(S));
f2 = @(w2) MDP_ObjectiveFunction(w2,std_all,S);
w0 = 1/p*ones(p,1);
weights = fmincon(f2,w0,[],[],UB,beq,LB',UB',[]);
end
end

