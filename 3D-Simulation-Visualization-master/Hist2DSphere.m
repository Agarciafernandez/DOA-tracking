function histY_n=Hist2DSphere(Y,varargin)
% HIST2DSPHERE Finds histogram counts on the HEALPix tesselated unit sphere. 
%              Data points (in standard cartesian format) are sorted, 
%              grouped and counted by which HEALPix pixel they are in.
%
% histY_n=Hist2DSphere(Y,nPix);
%%
% Examples of correct usage:
%
% nSide=2^2;
% nPix= nSide2nPix(nSide);
% kappa = 4.2;
% bet = 4.5 ;
% gamm = -3.5; 
% Psi= pi/2;
%
% Mu=[0 0 1];
% n = 2^12
%
% Y =Random_FB4(kappa,gamm,Mu,n);
%
%
% histY_n = Hist2DSphere(Y,nPix); 
% histY_n = Hist2DSphere(Y);
%%
%
% Required Inputs
% Y      n-by-3 dataframe, where each row is a datapoint in standard
%        cartesian coordinate format, and n is the sample size
%
% Optional Inputs
% nPix   number of HEALPix pixels (resolution parameter), should be of the 
%        form nPix = 3*4^k, for k = 1,2,...
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
dbl = {'double'};
 
% Check characteristics for the required and optional parameters
addRequired(p,'Y',@(x)validateattributes(x,dbl,{'ncols',3},'Y'));
addOptional(p,'nPix',[],@(x)validateattributes(x,dbl,{'scalar'},'nPix'));
 
p.parse(Y,varargin{:});
nPix = p.Results.nPix;

% Default Values
if isempty(nPix)
    nPix = 768;
end
if mod(log2(nPix/3),2)~= 0
    errormessage = 'Error:\nnPix must be of the form nPix=3*4^k \n for k=1,2,...';
    error('something:anything',errormessage)
end

%% Histogram Count Code
% normalized:  mean(histY)=1!!
nSide=nPix2nSide(nPix); % HEALPix tesselation parameter
Ypix=vec2pix(nSide,num2cell(Y,2)); % pixel numbers ,'nest',false 
Ysort=sort(Ypix); % Ysort(1:10) Ysort(end-10:end)
histY=zeros(nPix,1);
DYsort=diff(Ysort);
yd=find(DYsort~=0);
histY(Ysort(1)) = yd(1); % first one
ind1=cumsum(DYsort(DYsort~=0))+Ysort(1); % ind1(end-3:end)
Li=length(ind1)-1;
histY(ind1(1:Li))=diff(yd)'; % histY(1:3); histY(end-3:end)
histY(Ysort(end)) =length(Ysort)-yd(end); % last one
histY_n=histY*nPix/length(Y)/4/pi ; %;  integral of histogram is 1;