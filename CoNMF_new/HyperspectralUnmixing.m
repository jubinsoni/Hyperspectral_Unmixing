clear all ; close all; clc
%modified code required for transforming tif to mat file.
A= imread('clip_tif.tif');
[r,c,b]= size(A);
B= reshape(A, (r*c),b);
B=double(B');
save('clip_tif_matfile.mat','B');
l= load('clip_tif_matfile.mat');
x= size(l.B);
hdr = hdrread('reflectance1.hdr')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No_bands = 171;
% Lines= 316;
% Columns=103;
% load('clip_tif_matfile.mat'); % load the hyperspectral data
% 
% L = x(1);
% n=x(2);
% q = 14;   % estimated number of  endmembers
% 
% 
% verbose='on';
% noise_type = 'additive';        
% 
% [w Rn] = estNoise(B,noise_type,verbose);
% [kf,Ek,E]=hysime(B,w,Rn,1e-5,verbose);    % estimate the p
% 
% 
% Ek = E(:,1:q);
% Y = Ek'*B; UUp =Ek; 
% 
% 
% figure(1000)
% [Ahat,Xhat,rerr,L] =  ...
%        conmf(Y,q, ...
%             'POSITIVITY','yes', ...
%             'ALPHA',(1e-8)*sqrt(n), ...     %l_21 regularization parameter
%             'BETA',  10^(-3)*(n*q), ...      %minimum volume regularization parameter
%             'ADDONE','yes', ...
%             'AO_ITERS',30, ...
%             'DELTA', 1e-4,  ...    % (STOP) relative reconstruction error
%             'CSUNSAL_ITERS',100, ...
%             'MU_A', 1e-4*(n*q), ...  %proximity weight for A
%             'MU_X', 1e-1, ...      %proximity weight for X
%             'SPHERIZE', 'no', ...         %{'no','cov', 'M'}
%             'MIN_VOLUME', 'boundary', ...  %{'boundary', 'center', 'totalVar'} 
%             'VERBOSE','no', ...
%             'PLOT_A', 'yes');
% 
% 
%         
% figure(20)
% plot(UUp*Ahat)
% title( 'Ahat')
% 
% figure(40)
% plot(rerr)
% title('Reconstruction error')
% 
% 
% figure(50)
% plot(L)
% title('Objective function')
% 
% 
% 
% for i=1:q
%     figure(100+10*i)
%     imagesc(reshape(Xhat(i,:)',Lines,Columns));
% end
%     