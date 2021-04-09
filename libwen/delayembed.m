function X = delayembed(Y,tau,m)
% X = delayembed(Y, tau, m);
%
% Takens embedding with delay tau and embedding dimension m.
% Builds predictive coordinates so X=[y(t),y(t-tau),...,y(t-(m-1)*tau)]        
% 
% INPUT:
%   Y:   T x d dimensional observation matrix containing.
%   tau: Time delay to be used for embedding.
%   m:   Embedding dimension
%        
% OUTPUT:
%   X: (T-(m-1)*tau) x m state space matrix
%

    if nargin < 2
        error('Error: Needs at least three input variables!')
    end
    

%% Parameters for reconstruction
    T = length(Y);
    L = 1+(m-1)*tau;
  
    X = zeros((T-L+1),size(Y,2)*m);   
 
%% E-dimensional delay reconstruction with a delay of tau
    for t=1:(T-L+1)      
        dels = Y((L+t-1):-tau:(L+t-1-(m-1)*tau),:);
        X(t,:) = reshape(dels',1,size(Y,2)*m);
    end
end