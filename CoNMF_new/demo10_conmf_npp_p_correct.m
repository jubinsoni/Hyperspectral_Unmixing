%  demo10_npp_p_correct
%
%  This demo illustrates the collaborative hyperspectral unmixing via CoNMF
%  operating without  pure pixels and a correct subspace dimention (p = q).
%
%  ------------------------------------------------------------------------
%  SCENARIO:  non pure pixels,  correct  number of end, endmembers (p=q)
%   
%   3 < p <= 8
%   SNR >= 30 dB
%  ------------------------------------------------------------------------

%   
%        
%         
%  See
%
% J. Li, J. M. Bioucas-Dias, A. Plaza and L. Liu, "Robust Collaborative 
% Nonnegative Matrix Factorization for Hyperspectral Unmixing," 
% in IEEE Transactions on Geoscience and Remote Sensing, vol. 54, no. 10, 
% pp. 6076-6090, Oct. 2016. doi: 10.1109/TGRS.2016.2580702
%
%
%  -----------------------------------------------------------------------
%  
%  Collaborative Nonnegative Matrix Factorization (CoNMF):
%
%  Let  -> Y [L,np]  matrix containing the observed spectral vectors in its
%                   columns
%
%       -> M  [L,p] mixing matrix with p spectral signatures (usually L >> p)
% 
%       -> S  [p,np] abundance matrix  
%
%       -> N  {L,np] additive Gaussian Noise
%
%     
%  Linear mixing observation model     
%
%     Y = MS + N    with   S>= 0,   and   sum(S) = ones(1,np)
%
% -----------------------------------------------------------------------
%  Optimization problem: 
%
%   min  (1/2)  ||AX-Y||^2_F + alpha ||X||_{2,1} + beta ||BA-O||^2_F
%    X >= 0
%    A \in (q-1)affine set defined by the data
%
%   OPTIONAL CONSTRAINTS:
%
%       Sum-To-One sum(X) = ones(1,np)
%
%   where  ||X||_{2,1} = sum(sqrt(sum(X.^2,2))) is the l_{2,1} mixed 
%   norm which promotes rows of zeros in X and  ||BA-O||^2_F is a term
%   with minimum volume flavor.  The exact meaning of this term depends
%   in the couple(B,O):
%
%  'MIN_VOLUME' - {'boundary', 'center', }  define the minimum volume 
%                  term ||BA-O||^2_F
%
%                 'boundary'        B= I, O =  extremes given by VCA
%                 'center'          B=I,  O = center of mass
%                 'totalVar'        B = sqrtm(q*eye(q)-ones(q)/q/(q-1)), O = 0
%   
%
% Author: Jose Bioucas-Dias, February, November, 2014
%

clear all
close all

% define random states for reproducible result
%rand('state',23);
%randn('state',23);
%


%--------------------------------------------------------------------------
% Load Library (matrix A)
%--------------------------------------------------------------------------
% 1 - USGS           (L = 224; m = 498)
% 2 - USGS - pruned  (L = 224; m = 342) (3 deg)
% 3 - USGS - pruned  (L = 38;  m = 62) (10 deg)
% 4 - USGS - pruned  (L = 8;   m = 12) (20 deg)
% 5 - USGS - pruned  (L = 224; m = 60) (30 deg)
% 6 - iid Gaussian (0,1)
% 7 - iid Uniform [0,1]

% size of the mixing matrix [Lxm] (only applies to libraries 6 and 7) 
L = 200;

% set library
library = 3;    

% number of pixels
np = 4000;

% number of endmembers
p = 6;

% maximum number of endmembers per pixel
p_pix = 5;  % see description below

% SNR in dBs  (SNR = ||MX||^2_F / ||N||^2_F)    (Y = MS + N;)
SNR = 40

% abundance related parameters
MAX_PURIRY = 0.8;         % do not theshold abundances
OUTLIERS   = 0;           % no outliers
PURE_PIXELS = 'no';       % no pure pixels
SHAPE_PARAMETER = 1;      %(Dirichlet parameter) abundances uniformely 
                          %distribted  on each pidel

%% -------------------------  abundance generation ------------------  
% 
% Assumptions:
%
%   1)  the dataset contains p endmenbers
%   2)  number of endmembers per pixel is p_pix
%   3)  SHAPE_PARAMETER is the Dirichlet parameter applied to the sets of 
%       p_pix endmembers active at each pixel
%   4)  the groups of active pixels are sets  with p_mix components. The 
%       sets are mod([i,1+1,...,i+p_mix-1], p) for i=1,...,p
%    5) the groups have equal probability of being active
%

% Dirchlet parameters                          
pdf_pars = zeros(p, p);
aux = [ones(1,p_pix) zeros(SHAPE_PARAMETER,p-p_pix)];
for i=1:p
    pdf_pars(i,:) = circshift(aux',[i-1])';
end
% add weights
pdf_pars = [1/p*ones(p,1) pdf_pars];



%% --------------------------------------------------------------------------
%       Start simulation
%--------------------------------------------------------------------------

switch library
    case 1  % A_1
        load USGS_1995_Library.mat
        wavelengths = datalib(:,1);
        [dummy, indexes] = sort(wavelengths);
        A = datalib(indexes,4:end);
        names = names(4:end,:);
        clear datalib;
    case 2
        load USGS_pruned_3_deg.mat
        A = B;
        clear B;
    case 3
         load USGS_pruned_10_deg.mat
          %A = B(1:6:end,:);
          A = B;
         clear B;
    case 4
         load USGS_pruned_20_deg.mat
         %A = B(1:30:end,:);
         A = B;
         clear B;
    case 5
         load USGS_pruned_30_deg.mat
         A = B;
        clear B;
    case 6
        A = randn(L,2*L);
    case 7
        A = rand(L,2*L);
    otherwise
        disp('Unknown library')
end

[L,m] = size(A);   
% L = number of bands; m = total number of materials in the library
% normalize A
% nA = sqrt(sum(A.^2));
% A = A./repmat(nA,L,1);


%% 
% -------------------------------------------------------------------
%            Generate data 
% ------------------------------------------------------------------- 

% mixing matrix
index = randperm(m);
M = A(:,index(1:p));
   
% generate the data
n_aux  = 0;
Y = [];
Xaux = [];
N = [];
while n_aux < np  
[Ya,Xaaux,Na] = spectMixGen(M,np, ...
         'Source_pdf', 'Diri_mix', ...
         'pdf_pars',pdf_pars,...
         'max_purity',ones(1,p), ...
         'no_outliers',OUTLIERS, ...
         'pure_pixels', PURE_PIXELS, ...
         'violation_extremes',[1,1.2], ....
         'snr', SNR, ...
         'noise_shape','uniform');
        
     
  mask =   sum(Xaaux > MAX_PURIRY) ==  0;
  Y = [Y Ya(:,mask)];
  Xaux = [Xaux Xaaux(:,mask)];
  N = [N Na(:,mask)];
  n_aux = length(Y);
     
end 
     
Y = Y(:,1:np);
X = Xaux(:,1:np);
N = N(:,1:np);
   
     
    

%%
%--------------------------------------------------------------------------
% Optimization
%--------------------------------------------------------------------------

% q is the an overestimate of p
q = p;


%% scattergram centered at the mass center
figure(10)
ym = mean(Y,2);

Yc = Y-repmat(ym,1,np);
Cy = Yc*Yc'/np;
[U,S] = svd(Cy);

plot(U(:,1)'*Yc, U(:,2)'*Yc,'.');

Mc = M-repmat(ym,1,p);

hold on 
plot(U(:,1)'*Mc(:,1:p), U(:,2)'*Mc(:,1:p),'ro','LineWidth',3);
plot(0, 0,'ko','LineWidth',3);
axis  on
hold off 
title('Projection of Y onto the firt two PCA components')

figure(1000)
[Ahat,Xhat,rerr,L] = ...
        conmf(Y,q, ...
            'POSITIVITY','yes', ...
            'ALPHA', 10^(-9)*sqrt(np/p), ...    % l_21 regularization parameter
            'BETA', 10^(-SNR/10-1)*(np*p) ,...    % minimum volume regularization parameter
            'ADDONE','yes', ...
            'AO_ITERS', 100, ...
            'DELTA', 1e-8,  ...            % (STOP) relative reconstruction error
            'CSUNSAL_ITERS',100, ...
            'MU_A', 1e-4*(np*p), ...       %proximity weight for A
            'MU_X', 1e-1, ...              %proximity weight for X
            'SPHERIZE', 'cov', ...         %{'no','cov', 'M'}
            'MIN_VOLUME', 'totalVar', ...  %{'boundary', 'center', 'totalVar'} 
            'TRUE_M',M, ... 
            'VERBOSE','no', ...
            'PLOT_A', 'yes');

figure(30)
imagesc(Xhat)
title('Estimated abundances')


% find alignment
P = align_matrices(Ahat,M,'angle');

figure(35)
plot([Ahat*P - M])
title('Error (Ahat-M)')

% compute SAD
aux = sum((Ahat*P).*M)./ sqrt(sum((Ahat*P).^2).*sum(M.^2));
fprintf('\n\nSAD (deg) = %2.2f\n', mean(abs(acos(aux)))*180/2/pi);


fprintf('\nRERR = %2.2f\n', rerr(end))


figure(45)
plot(L)
title('Objective function')

figure(50)
plot(rerr)
title('Reconstruction error')

error =sum((Y-Ahat*Xhat).^2);
figure(60)
hist(error,100)
title('histpgram of ||Y-AX||^2')


figure(65)
plot([M Ahat])
title('M and Ahat')


fprintf('\nNormalized RERR = %f\n', rerr(end)/norm(Y,'fro'))

