function h=Plot_DataRandomS2(Y)
% PLOT_DATARANDOMS2 Plots data points on the unit sphere (S2)
%
% h=Plot_DataRandomS2(Y)
%%
% Examples of correct usage:
%
% n=2^12; 
% Y = Random_Uni_Inv(n);
% h=Plot_DataRandomS2(Y);
%%
%
% Required Inputs
% Y     n-by-3 matrix/dataframe, where each row is a datapoint in standard 
%       cartesian coordinate format, and N is the sample size
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
dbl = {'double'};

% Check characteristics for the required and optional parameters
addRequired(p,'Y',@(x)validateattributes(x,dbl,{'ncols',3},'Y'));

%%
Rp=Y'; % convert form to (1:3,:) 
% figure    
resolution = 50; % graphical parameter
c = ones(resolution); 
[x,y,z] = sphere(resolution);
h2 = surf(x,y,z,1*c);
set(h2,'LineStyle','none');
colormap([0 1 0])
axis tight; hold on
% plot3(x1,y1,z1,'.','MarkerSize',15,'color',[0,0,0]) %nonuniform
h=plot3(Rp(1,:)',Rp(2,:)',Rp(3,:)','.','MarkerSize',4,'color',[0,0,0]); %nonuniform
hold off
view(3);
alpha(h2,.9);