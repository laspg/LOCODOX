% nlinfitci fits the model specified by MODELFUN to variables in the
% dataset
%
% SYNTAX
% [beta,mdl] = nlinfitci(X,Y,optfun,beta0,options)
%
% DESCRIPTION
% Fits the model specified by MODELFUN to variables in the dataset/table
% array DS, and returns the nonlinear model NLM. The coefficients are
% estimated using an iterative procedure starting from the initial values
% in the vector BETA0.
%
% INPUT
%   X         A design matrix of predictor (independent variable)
%             values, with one row for each value in Y and one column for
%             each coefficient.  However, X may be any array that MODELFUN
%             is prepared to accepted.
%   Y         A vector of response (dependent variable) values
%   optfun    A function, can be either of the following:
%                 1. A function, specified using @, that accepts two
%                 arguments, a coefficient vector and an array of predictor
%                 values, and computes a vector of Y values. 
%                 2. A text string, such as 'y ~ b0+b1*sin(b2*X)', defining
%                 the response and a mathematical expression for the
%                 function. Symbols in the function expression that match
%                 variable names in the dataset/table are taken to be
%                 predictors, and the others are coefficients to be
%                 estimated.
%   beta0     The initial values of the coefficient of the model
%
%  OUTPUT
%   beta      The confidence intervals for the parameter in the nonlinear
%             regression
%   mdl       The nonlinear model
%
% CALL :
%
% SEE ALSO
%   

% HISTORY
% according to recommendations by SCOR WG 142 "Quality Control Procedures
% for Oxygen and Other Biogeochemical Sensors on Floats and Gliders"
%
% Henry Bittig
% Laboratoire d'Oceanographie de Villefranche-sur-Mer, France
% bittig@obs-vlfr.fr
% 28.10.2015
%   $created: 28/10/2015 $author: Henry Bittig, 
%               Laboratoire d'Oceanographie de Villefranche-sur-Mer, France
%               bittig@obs-vlfr.fr
%   $Revision: version $Date: $author:
%   v1.2 25/05/2016   Emilie Brion, ALTRAN OUEST
%                     Use fitnlm instead of nlinfit. Output the model, not
%                     only the coefficients.
%   v1.3 04/07/2016   Emilie Brion,desc ALTRAN OUEST
%                     Update the descriptive header



function [beta,mdl] = nlinfitci(X,Y,optfun,beta0)

% Coefficients initialisation
beta0 = beta0(:);

% Get the model
mdl = fitnlm(X,Y,optfun,beta0);
beta(:,1)= mdl.Coefficients.Estimate;
% Confidence intervals for parameters in nonlinear regression
ci = nlparci(mdl.Coefficients.Estimate,mdl.Residuals.Raw,'covar',mdl.CoefficientCovariance);
beta(:,2) = (ci(:,2)-ci(:,1))/2;

% OLD version
% if nargin<5
%     options=[];
% end
% [beta,r,J,COVB,mse]=nlinfit(X,Y,optfun,beta0,options);
% % Confidence intervals for parameters in nonlinear regression
% ci=nlparci(beta,r,'covar',COVB);
% beta(:,2)=(ci(:,2)-ci(:,1))/2;
