function [A,X,rerr,L] = conmf(Y,q,varargin)

%% [A,X,rerr,L] = conmf(Y,q,varargin)
%
%  Collaborative nonnegative factorization (CoMNF)
%
%  
%
%% --------------- Description ---------------------------------------------
%
%  Collaborative Nonnegative Matrix Factorization (CoNMF):
%
%  Let  -> Y [L,np]  matrix containing the observed spectral vectors in its
%                    columns
%
%       -> M  [L,q] mixing matrix with p spectral signatures (usually L >> q)
%
%       -> S  [q,np] abundance matrix
%
%       -> N  {L,np] additive Gaussian Noise
%
%
%  Linear mixing model
%
%     Y = MS + N    with   S>= 0,   and   sum(S) = ones(1,np)
%
%  CoNMF solve the follwoing ptimization problem:
%
%   min  (1/2)  ||AX-Y||^2_F + alpha ||X||_{2,1} + beta ||BA-P||^2_F
%    X >= 0
%    A \in (q-1) affine set defined by the data
%
%    OPTIONAL CONSTRAINTS:
%
%       Sum-To-One sum(X) = ones(1,np)
%
%   where  ||X||_{2,1} = sum(sqrt(sum(X.^2,2))) is the l_{2,1} mixed
%   norm with promotes rows of zeros in X, and  the couple (B,O) 
%   define a minimum volume-type  regularizer:
%   
%   'MIN_VOLUME' - {'boundary', 'center', 'total variance'}
%   
%      'boundary' - > B=I, O = (VCA solution)
%    
%      'center' - >   B=I, O = (sample mean value)
%
%      'totalVar'     B = sqrtm(q*eye(q)-ones(q)/q/(q-1)), O = 0
%           
%       In the the three above cases,  the term beta ||BA-O||^2_F
%       pushes the columns of A  towards the simplex thus having a minimum 
%       volume flavor
%
%
%% -------------------- Line of Attack  -----------------------------------
%
%  CoNMF implemets proximal alternating sptimization  (PAO). See [1]
%
%  The 
%
% ------------------------------------------------------------------------
%%  ===== Required inputs =============
%
% Y - matrix with  L(channels) x N(pixels).
%     each pixel is a linear mixture of p endmembers
%     signatures y = M*x + noise,
%
%     CoNMF assumes that y belongs to an affine space. It may happen,
%     however, that the data supplied by the user is not in an affine
%     set. For this reason, the first step this code implements
%     is the estimation of the affine set the best represent
%     (in the l2 sense) the data.
%
% q - overestimate of the number of independent columns of M. Therefore,
%   M spans a (q-1)-dimensional affine set.
%
%
%%  ====================== Optional inputs =============================
%
%  'POSITIVITY'  = {'yes', 'no'}; Enforces the positivity constraint:
%                   X >= 0
%                   Default 'yes'
%
%  'ADDONE'  = {'yes', 'no'}; Enforces the positivity constraint: X >= 0
%              Default 'no'
%
%
%  'ALPHA' - regularization parameter for the l_{2,1} norm
%
%  'BETA'  - regularization parameter for the minimum volume term
%
%  'MU_A' -  weight for A  proximal term
%
%  'MU_X' -  weight for X  proximal term
%
%  'AO_ITERS' - CoNMF  number of iterations
%
%  'DELTA' -  relative reconstruction error. Used as STOP criterion
%
%  'THETA' -  thresold to detect the active rows of X
%             default = 1e-3*n;
%
%  'CSUNSAL_ITERS' - CSUNSAL number of iterations
%                    Default: 100;
%
%  'SPHERIZE' -  {'no','cov', 'M'} sperize the observed data Y
%                'cov'  based on the covariance matrix
%                'M'    based on an estimated mixing matrix
%   
%  'MIN_VOLUME' - {'boundary', 'center', }  define the minimum volume 
%                  term ||BA-O||^2_F
%
%                 'boundary'    B= I, O =  extremes given by VCA
%                 'center'      B=I,  O = center of mass
%                 'totalVar'    B = sqrtm(q*eye(q)-ones(q)/q/(q-1)), O = 0
%
%  'A0' -  initial A
%
%  'TRUE_M' - accept the true mixing matrix (only for displaying purposes)
%
%  'VERBOSE'   = {'yes', 'no'};
%                 'no' - work silently
%                 'yes' - display warnings and plot figures
%                  Default 'no'
%
%
%
%
%%  =========================== Outputs ==================================
%
% A  =  [Lxq] estimated mixing matrix
%
% X  =  [qxnp] estimated abundances
%
% rerr ->  reconstruction error along the iterations
%
% L   ->   objective fuction along the iterations
%
%% -------------------------------------------------------------------------
%
% Copyright (November, 2014):   Jose Bioucas-Dias (bioucas@lx.it.pt)
%                               Jun li (lijun48@mail.sysu.edu.cn)
%
% CoNMF is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
%
% -------------------------------------------------------------------------

%%
%--------------------------------------------------------------
% test for number of required parametres
%--------------------------------------------------------------
if (nargin-length(varargin)) ~= 2
    error('Wrong number of required parameters');
end
% data set size
[d,n] = size(Y);
if (d<q)
    error('Insufficient number of columns in Y');
end
%%
%--------------------------------------------------------------
% Set the default for the optional input parameters
%--------------------------------------------------------------
%positivity
positivity = 'yes';
%addone
addone = 'no';
%regularization parameters
lam_x = 1e-6*sqrt(n/1000); % mixed norm  regularization parameter
lam_a = 5e-1*(n/1000);      % minumum volume regularization parameter
% parameter for the quadratic proximal term for A
mu_a = 1e0*(n/1000);
% parameter for the quadratic proximal term for X
mu_x = 1e0;
% maximum number of AO iterations
AOiters = 200;
% relative reconstruction error (RECONS_ERR) for the termination test
delta = 1e-4;
% thresold to detect the active rows of X
theta = 1e-3*n;
%  CSUNSAL number of iterations
csunsal_iters = 100;
% true A
M_true = [];
% display only CoNMF warnings
verbose = 'no';
%spherize
spherize = 'M';
% minimum volume tetrm
min_volume = 'boundary'
% no initial simplex
A = [];
% use csunsal_hinge instead of csunsal
hinge = 0;
% plot evolution of A(t)
plot_a_t = 'no';

%--------------------------------------------------------------
% Read the optional input parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional input parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'POSITIVITY'
                positivity = varargin{i+1};
            case 'HINGE'
                hinge = varargin{i+1};
            case 'ADDONE'
                addone = varargin{i+1};
            case 'ALPHA'
                lam_x = varargin{i+1};
            case 'BETA'
                lam_a = varargin{i+1};
            case 'MU_A'
                mu_a = varargin{i+1};
            case 'MU_X'
                mu_x = varargin{i+1};
            case 'AO_ITERS'
                AOiters = varargin{i+1};
            case 'DELTA'
                delta = varargin{i+1};
            case 'THETA'
                theta = varargin{i+1};
            case 'CSUNSAL_ITERS'
                csunsal_iters = varargin{i+1};
            case 'TRUE_M'
                M_true = varargin{i+1};
            case 'SPHERIZE'
                spherize = varargin{i+1};
            case 'MIN_VOLUME'
                min_volume = varargin{i+1};
            case 'A0'
                A = varargin{i+1};
            case 'VERBOSE'
                verbose = varargin{i+1};            
            case 'PLOT_A'
                 plot_a_t = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end;
    end;
end

%%
%--------------------------------------------------------------
% take care of output variable
%--------------------------------------------------------------

% compute objective function
comp_obj_func = 0;

if nargout == 4
    comp_obj_func = 1;
end

%%
%--------------------------------------------------------------
% Local variables
%--------------------------------------------------------------

% regulatization parameter for spherization
lam_sphe_reg  = 1e-6;
% VCA number of VCA runs
VCA_runs = 50;

%%
%--------------------------------------------------------------------
% identify the affine space that best represent the data set y in R^q
%--------------------------------------------------------------------
my = mean(Y,2);
Y = Y-repmat(my,1,n);
[Uq,D] = svd(Y*Y'/n);
Uq = Uq(:,1:q-1);    % I'm using svd because I've found errors in svds
sing_values = diag(D);

% define affine projection of dimension q-1 and lift
aff_proj = @(Z) [ Uq'*Z; ones(1,size(Z,2))];

% compute the coordinates of Y (centered) wrt Uq and lift
Y = aff_proj(Y);

% ------------------------------------------------------------------------
% from here on  we are are representing  the data as
%
%   y = [yp ; 1];
%
% where yp = Uq'(y-my) \in R^(q-1)
%
% -----------------------------------------------------------------------


%%
%------------------------------------------
% spherize if requested
%------------------------------------------
if strcmp(spherize,'cov')
    C = diag([1./(sqrt(sing_values(1:q-1))+lam_sphe_reg); 1]);
elseif strcmp(spherize,'M')
    % estimate M with VCA
    Aux = VCA(Y,'Endmembers',q,'SNR',1,'verbose',verbose);
    C = inv(Aux + lam_sphe_reg*eye(q));
else
    C = eye(q);
end

% spherize data
Y = C*Y;
% represent A in the new coordinate system
if ~isempty(A)
   A = C*aff_proj(A-my*ones(1,q));
end
% represent M_true in the new coordinate system
if ~ isempty(M_true)
   M_true = C*aff_proj(M_true-my*ones(1,size(M_true,2)));
end


%--------------------------------------------------------------------
% NOTE: the coordinates the affine set defined by Y may  be written as 
%
%            y = yc + F*w
%
% where     yc = C*[0, ...., 0, 1]'  and F = orth(C*[eye(q-1);ones(1,q-1)])
%% 
%--------------------------------------------------------------
yc = C*[zeros(q-1,1); 1];
F = [eye(q-1); zeros(1,q-1)];
% orthogonalize
F = orth(C*F);


%%
%---------------------------------------------
%  Initialization
%---------------------------------------------
small = 1e-4;
vol_ant = -inf;

% Initialize with VCA
for i=1:VCA_runs
    Aux = VCA(Y,'Endmembers',q,'SNR',1,'verbose',verbose);
    %Aux = svmax(Y,'Endmembers',q,'SNR',1,'verbose',verbose);
    vol = sum(log(svd(Aux) + small ));
    if vol > vol_ant
        Avca = Aux;
        vol_ant = vol;
    end
end
if isempty(A)
    A = Avca;
    % expand Q
    Am = mean(A,2);
    Am = repmat(Am,1,q);
    dQ = A - Am;
    % fraction: multiply by p is to make sure Q0 starts with a feasible
    % initial value.
    A = A + 0.0*q*dQ;
end



% initial abundances
X = sunsal(A,Y,'POSITIVITY',positivity,'VERBOSE',verbose,'ADDONE',addone, ...
    'lambda', 1e-3,'AL_ITERS',1000, 'TOL', 1e-8, 'verbose',verbose);

% A is represented in the affine set as A = Yc + F*W
Yc = repmat(yc,1,q);

% set the couple(B,O) for the minimum volume term
if strcmp(min_volume, 'boundary')
    %use VCA extremes  
    O = Avca/q;
    B = eye(q)/q;
elseif strcmp(min_volume, 'center')
    % mass-center
    O = repmat(yc,1,q)/q;
    B = eye(q)/q;
else  % total variance
    O = zeros(q);
    B= ((q*eye(q)-ones(q))/q/(q-1));
    [Uaux,Saux ] = svd(B);
    B=Uaux*sqrt(Saux)*Uaux';
end

% use vca or mass-center as points close to vertices to define the 
% minimum volume term

 
if strcmp(plot_a_t, 'yes')
    
    plot(Y(1,:), Y(2,:),'.');
    hold on
    plot(Avca(1,:), Avca(2,:),'go','LineWidth',3);
    plot(yc(1),yc(2),'ko','LineWidth',3);
    
    if ~ isempty(M_true)
        plot(M_true(1,:),M_true(2,:)','ro','LineWidth',3);
    end
end


%---------------------------------------------
%   Constants
%---------------------------------------------

K = lam_a*F'*B'*B*F + mu_a*eye(q-1);
[Uk,Sk] = svd(K);

Sk = repmat(diag(Sk),1,q); 

Q = lam_a*B'*(B*Yc-O);


% inilialize rerr and L
% norm of Y
norm_y = norm(Y,'fro');
% initial reconstruction errors
rerr(1) = norm(Y-A*X,'fro');

L(1)=0.5*norm(Y-A*X,'fro')^2 + lam_x*sum(sqrt(sum(X.^2,2))) + ...
            0.5*lam_a*norm(B*A-O,'fro')^2;
        
        
%---------------------------------------------
%   main body
%---------------------------------------------       
t=1;
k=1;

while (t <= AOiters) && (rerr(t)/norm_y > delta)

    if strcmp(verbose, 'no')
        fprintf('OA iter = %d \n',t);
    end
    
    XXT = X*X';
    % update A
    [Uh,Sh] = svd(XXT);    
    G = -F'*(Yc*XXT - Y*X' + Q + mu_a*(Yc -A));
    Dkh = 1./(repmat(diag(Sh)',q-1,1) + Sk);
    W = Uk*(((Uk'*G*Uh).*Dkh))*Uh';
   
    % compute A
    A = Yc + F*W;
    
    rerr(t+1) = norm(Y-A*X,'fro');

    if comp_obj_func
        L(t) = 0.5*norm(Y-A*X,'fro')^2 + lam_x*sum(sqrt(sum(X.^2,2))) + ...
            0.5*lam_a*norm(B*A-O,'fro')^2;
    end
    
    % update X
    Y_aux = [Y; sqrt(mu_x)*X];
    A_aux = [A; sqrt(mu_x)*eye(q)];
    
    X = clsunsal(A_aux,Y_aux,'POSITIVITY',positivity,'VERBOSE',verbose,'ADDONE',addone, ...
        'lambda', lam_x,'AL_ITERS',csunsal_iters, 'TOL', 1e-6, ...
        'X0',X,'verbose', verbose);
    
%     X = clsunsal_hinge(A_aux,Y_aux,'HINGE',hinge,'VERBOSE',verbose,'ADDONE',addone, ...
%         'lambda', lam_x,'AL_ITERS',csunsal_iters, 'TOL', 1e-6, ...
%         'X0',X,'verbose', verbose);
    
 
    
    
    if strcmp(plot_a_t, 'yes')
        plot(A(1,:), A(2,:),'k.','LineWidth',3);
        drawnow;
    end
    t = t+1;
    
end % main body


if strcmp(plot_a_t, 'yes')
    plot(A(1,:), A(2,:),'ms','LineWidth',2);
end

hold off
% undo spherization and back to the original representation
A = inv(C)*A;
A = Uq*A(1:q-1,:);
A = A + repmat(my,1,q);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
