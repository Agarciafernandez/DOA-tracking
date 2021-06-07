function Y = Random_Uni_Norm(varargin)
% RANDOM_UNI_NORM Generates a random sample of size n from the Uniform 
% distribution in S2
%
% Y = Random_Uni_Norm(n)
%%
% Examples of correct usage:
%
% n=2^12; 
%
% Examples of correct function construction: 
% Y = Random_Uni_Norm(n)
% Y = Random_Uni_Norm
%%
%
% Optional Inputs
% n           sample size (where n is a positive integer)
%
% Authors: Gy.Terdik, B.Wainwright
%
% Copyright 2018 Gy.Terdik
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
%% Verify/Default parameter values

% Required parameter characteristics
p = inputParser;
class = {'numeric'};
pos_int = {'integer','positive'};

% Check characteristics for the required and optional parameters
addOptional(p,'n',[],@(x)validateattributes(x,class,pos_int));

p.parse(varargin{:});
n = p.Results.n;

% Default Values
if isempty(n)
    n = 2^15;
end

%% density
R1 = normrnd(0,1,3,n); % X1_i,X2_i,X3_i iid Normal(0,1) for i=1,2,...,n
R1norm = sqrt(sum(R1.^2,1));
R1norm = kron(ones(3,1),R1norm);
Y = R1'./R1norm'; % normalize vectors, so all are of length 1
return