function [Zw, whiteningMatrix, dewhiteningMatrix] = whitenRows(Z,maxEig)
%
% Syntax:       [Zw, T] = whitenRows(Z,maxEig);
%               
% Inputs:       Z is an (d x n) matrix containing n samples of a
%               d-dimensional random vector
%
%               maxEig is a scalar determining the number of singular
%               vectors to be used.
%               
% Outputs:      Zw is the whitened version of Z
%               
%               T is the (d x d) whitening transformation of Z
%               
% Description:  Returns the whitened (identity covariance) version of the
%               input data
%               
% Notes:        (a) Must have n >= d to fully whitenRows Z
%               
%               (b) Z = T \ Zcw
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         November 1, 2016
%

% covarianceMatrix = cov(Z', 1);
% [E, D] = svd (covarianceMatrix);
% %reduce dime here
% %whiten
% whiteningMatrix = inv (sqrt (D)) * E';
% dewhiteningMatrix = E * sqrt (D);
% newVectors =  whiteningMatrix * vectors;

% Compute sample covariance
R = cov(Z');

% Whiten data
[U, S, ~] = svd(R,'econ');

%added by Tosif
S=S(1:maxEig,1:maxEig);U=U(:,1:maxEig);
whiteningMatrix  = diag(1 ./ sqrt(diag(S))) * U';
dewhiteningMatrix = U * sqrt (S);
Zw = whiteningMatrix * Z;

%commented out the following
% T  = U * diag(1 ./ sqrt(diag(S))) * U';
% Zw = T * Z;
