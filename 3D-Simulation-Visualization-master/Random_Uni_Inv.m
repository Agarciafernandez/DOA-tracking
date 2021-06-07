function Y = Random_Uni_Inv(varargin)
% RANDOM_UNI_INV Generates a random sample of size n from the Uniform 
%                distribution in S2
%
% Y = Random_Uni_Inv(n)
%%
% Examples of correct usage:
%
% n=2^12; 
%
% Examples of correct function construction: 
% Y = Random_Uni_Inv(n)
% Y = Random_Uni_Inv
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
U = rand(n,1);
cX = 2*U-1; % cos(theta)
psi = 2*pi*rand(n,1);
sX = sqrt(1-cX.^2); % sin(theta)
% polar coordinates: 
Y = [cos(psi).*sX,sin(psi).*sX,cX];
return