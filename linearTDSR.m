function results = linearTDSR(X,r,u,alpha,opto,optotype)
    
    % Linear TD-SR algorithm
    %
    % USAGE: results = linearTDSR(X,r,opto)
    %
    % INPUTS:
    %   X - [N x D] stimulus sequence (N = # trials, D = # features)
    %   r - [N x 1] reward sequence
    %   u (optional) - initial values of the features (default: all zeros)
    %   alpha (optional) - 2 element vector with the SR Matrix alpha first,
    %                      and the TDRL alpha second (default: 0.06, 0.03)
    %   opto (optional) - sequence of optogenetic perturbations (default: all zeros)
    %   optotype (optional) - string of 'none' (default) or 'tonic' for
    %                         prolonged stimulation or 'phasic' for short
    %                         stimulation
    %
    % OUTPUTS:
    %   results - [N x 1] structure with the following fields:
    %               .R - reward estimate
    %               .V - value estimate
    %               .dt - prediction error
    %               .W - weight matrix
    %
    % Sam Gershman, Dec 2017, edited by Matt Gardner April 2018
    
    [N,D] = size(X);
    X = [X; zeros(1,D)];% add buffer at the end
    W = zeros(D);
    r = [r; 0];
    
    if nargin < 3
        u = zeros(D,1); alpha_W = 0.2; alpha_u = 0.1; opto = zeros(N,1); optotype = 'none';
    elseif nargin < 4
        u = zeros(D,1); alpha_W = alpha(1); alpha_u = alpha(2); opto = zeros(N,1); optotype = 'none';
    elseif nargin < 5
        alpha_W = alpha(1); alpha_u = alpha(2); opto = zeros(N,1); optotype = 'none';
    elseif nargin < 6
        disp('the optotype input is missing')
        return
    else
        alpha_W = alpha(1); alpha_u = alpha(2); 
    end
    
    gamma = 0.95;
    
    %For these experiments, we're only using currently active elements for
    %opto stimulation such that lambda = 0
    lambda = 0;
    E = 0;
    
    for n = 1:N
         
        switch optotype
            
            case 'none'
                dt = X(n,:) + (gamma*X(n+1,:) - X(n,:))*W;
            
            case 'tonic'
                dt = (1 + opto(n))*(X(n,:) + (gamma*X(n+1,:) - X(n,:))*W);
       
            case 'phasic'
                E = lambda*E + X(n,:);
                dt = X(n,:) + (gamma*X(n+1,:) - X(n,:))*W + opto(n)*E;  
          
        end
        
        %This gets the error for U
        du = r(n) - X(n,:)*u; 
        
        %This corrects for the excitatory/inhibitory assymetry by reducing the inhibitory components
        %Inh to Exc ratio currently set at 1:4
        dt(dt < 0) = dt(dt < 0)/4;
        du(du < 0) = du(du < 0)/4;        
        
        % store results
        results(n).R = X(n,:)*u;
        results(n).V = X(n,:)*W*u;
        results(n).dt = dt;
        results(n).W = W;
        results(n).U = u;
        results(n).AllV = W*u; 
        
        % update
        W = W + alpha_W*X(n,:)'*dt;
        u = u + alpha_u*X(n,:)'*du;
        
    end