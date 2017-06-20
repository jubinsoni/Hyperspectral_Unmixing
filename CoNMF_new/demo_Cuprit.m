%  demo_hyperspectral
%
%  This demo illustrates the collaborative hyperspectral unmixing via CoNMF
%  for the hyperspectral data.
%        
%         
%  
%  -----------------------------------------------------------------------
%

clear all; close all; clc

No_bands = 188;
load cuprite_ref.mat % load the cuprite data

[L,n] = size(x);

q = 11;   % estimated number of  endmembers


verbose='on';
noise_type = 'additive';        

[w Rn] = estNoise(x,noise_type,verbose);
[kf,Ek,E]=hysime(x,w,Rn,1e-5,verbose);    % estimate the p


Ek = E(:,1:q);
Y = Ek'*x; UUp =Ek; 


figure(1000)
[Ahat,Xhat,rerr,L] =  ...
       conmf(Y,q, ...
            'POSITIVITY','yes', ...
            'ALPHA',(1e-8)*sqrt(n), ...     %l_21 regularization parameter
            'BETA',  10^(-3)*(n*q), ...      %minimum volume regularization parameter
            'ADDONE','yes', ...
            'AO_ITERS',30, ...
            'DELTA', 1e-4,  ...    % (STOP) relative reconstruction error
            'CSUNSAL_ITERS',100, ...
            'MU_A', 1e-4*(n*q), ...  %proximity weight for A
            'MU_X', 1e-1, ...      %proximity weight for X
            'SPHERIZE', 'no', ...         %{'no','cov', 'M'}
            'MIN_VOLUME', 'boundary', ...  %{'boundary', 'center', 'totalVar'} 
            'VERBOSE','no', ...
            'PLOT_A', 'yes');


figure(20)
plot(UUp*Ahat)
title( 'Ahat')

figure(40)
plot(rerr)
title('Reconstruction error')


figure(50)
plot(L)
title('Objective function')



for i=1:q
    figure(100+10*i)
    imagesc(reshape(Xhat(i,:)',Lines,Columns));
end
    