% DOXY_corr_linreg computes the linear regression between variable and the
% reference.
%
% SYNTAX
%    [LRfit, LRFullFit, LRcoef] = DOXY_corr_linreg(myVar,myVarRef,isok,linReg,fitted0,okReg)
%
% DESCRIPTION
% DOXY_corr_linreg computes the linear regression between variable and the
% reference. adapts according to the availability of toolbox matlab
% (statistics or curve fitting toolbox).
%
% INPUT
%
%  myVar    (array)      variable field
%  myVarRef (array)      variable reference field
%
%  isok (logical array)  value in myVar an myVarRef to use to compute the
%                        linear regression
%
%  linReg (struct)       regression equation description, and Z intercept
%                        information. linReg.intercept = false means that
%                        the regression equation is forced to origin.
%                        Example: linReg = {'y ~ a*x + b','y ~ a*x'};
%                         linReg =
%                             equation: {'y ~ a*x + b'  'y ~ a*x'}
%                             intercept: [1 0]
%
%  fitted0 (logical)     (Optional) use the minimum in data (0) or the
%                        value O (1) to build the fitted data. Default
%                        = 0 (use the minimum in data).
%
%  okReg (double)        (Optional) Index in linReg to select the regression
%                        equation to use. by default, get all the equation
%                        available in linReg.
%
% OUTPUT
%   LRfit (array)         array of fitted data (min and max) : Line.
%                         Size = [length(okReg),2]
%
%   LRFullFit (array)     array of fully fitted data : all data fitted with
%                         the linear regression coefficient.
%                         Size = [length(okReg),size(myVar,1),size(myVar,2)]
%
%   LRcoef (array)        array of linear regression coefficients and
%                         R-squared indicator.
%                         Size = [length(okReg),3];
%
% CALL :
%
% SEE ALSO
%   DOXY_corr_compute, DOXY_corr_main
%
% HISTORY
%   $created: 02/07/2018 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:
%   v1.1 05/07/2018   Emilie Brion, Altran Ouest
%                     myVar(isok) and myVarRef(isok) should be LINE array
%                     sse correctd : no more b*data, but a*data
%                     warndlg : mode change to "replace" instead of
%                     "unmodal".
%                     license for toolbox : the check is reliable.
%   05/03/2019        Marine Gallian, Altran Ouest
%                     Fix a bug appearing when the statistic toolbox is
%                     not available

function [LRfit, LRFullFit, LRCoef] = DOXY_corr_linreg(myVar,myVarRef,isok,linReg,fitted0,okReg)


% =========================================================================
%% Initialisation
% =========================================================================
if nargin == 4
    okReg = 1:length(linReg.intercept);
    fitted0 = 1;
end
if nargin  == 5
    okReg = 1:length(linReg.intercept);
end

if fitted0 == 0
    minFit = nanmin(nanmin(myVar));
else
    minFit = 0;
end

[status,errmsg] = license('checkout','statistics_toolbox');
if status ~= 1
    [status,errmsg] = license('checkout','curve_fitting_toolbox');
end
status = license('inuse');
status = {status.feature};

% =========================================================================
%% COmpute Linear Regression, depending on matlab toolbox availability
% =========================================================================
warning('off','stats:LinearModel:RankDefDesignMat'); %Supress Warning indicating that there is several regression possible for the dataset
for i = 1:length(okReg)
    if any(~cellfun('isempty',strfind(status,'statistics')))
        if linReg.intercept(okReg(i)) == true
            mdl = fitlm(myVar(isok),myVarRef(isok),'linear','Intercept',true);
            a = mdl.Coefficients.Estimate(2);
            b = mdl.Coefficients.Estimate(1);
        else
            mdl = fitlm(myVar(isok),myVarRef(isok),'linear','Intercept',false);
            a = mdl.Coefficients.Estimate;
            b = 0;
        end
        fitted = polyval([a b],[minFit nanmax(nanmax(myVar))]);
        LRfit(i,:) = fitted;
        LRFullFit(i,:,:) = polyval([a b],myVar);
        LRCoef(i,:) = [a b mdl.Rsquared.Ordinary];
    elseif  any(~cellfun('isempty',strfind(status,'curve_fitting_toolbox')))
        equation = linReg.equation{okReg(i)};
        equation = regexp(equation,'y ~ ','split');
        equation = equation{2};
        % X must be a matrix with one or two columns.
        if ~iscolumn(myVar(isok))
            myVarOk = myVar(isok)';
        else
            myVarOk = myVar(isok);
        end
        if ~iscolumn(myVarRef(isok))
            myVarOkRef = myVarRef(isok)';
        else
            myVarOkRef = myVarRef(isok);
        end
        [f,gof,~] = fit(myVarOk,myVarOkRef,equation);
        regCoeff = coeffvalues(f);
        if length(coeffvalues(f)) == 1
            regCoeff(end+1) = 0;
        end
        a = regCoeff(1);
        b = regCoeff(2);
        fitted = polyval([a b],[minFit nanmax(nanmax(myVar))]);
        LRfit(i,:) = fitted;
        LRFullFit(i,:,:) = polyval([a b],myVar);
        LRCoef(i,:) = [a b gof.rsquare];
    else
        warndlg({'Curve Fitting and Statistics Toolbox Licenses unavailable : Matlab can''t apply your regression equation with native function.';... The correction will be compute with a polynomial equation'}]
            '';...
            'The correction will be simply compute with the polynomial equation y ~ a*x, and rsquare is simply computed.';''},...
            'DOXY_corr_compute : Toolbox unavailable','replace')
        if linReg.intercept(okReg(i))== true
            [a,b,rsquare,~,~] = lsqfitma(myVar(isok),myVarRef(isok));
        else
            % compute the slope (y ~ a*x) : the array should be line
            % array, and the second dimension should be the same
            refSize = size(myVarRef(isok));
            varSize =  size(myVar(isok));
            if refSize(1) ~= 1
                lineVarRef = myVarRef(isok)';
            else
                lineVarRef = myVarRef(isok);
            end
            if varSize(1) ~= 1
                lineVar = myVar(isok)';
            else
                lineVar = myVar(isok);
            end
            a = lineVarRef/lineVar;
            b = 0;
            data = myVar(isok); refdata = myVarRef(isok);
            % Determine the size of the vector
            n = length(data);
            % Calculate sums and other re-used expressions
            Sdata = sum(data);
            databar = Sdata/n;
            sse = sum((a*data - refdata).^2);
            sst = sum((data - databar).^2 );
            rsquare = 1 - sse/sst;
        end
        fitted = polyval([a b],[minFit nanmax(nanmax(myVar))]);
        LRfit(i,:) = fitted;
        LRFullFit(i,:,:) = polyval([a b],myVar);
        LRCoef(i,:) = [a b rsquare];
    end
end
warning('on','stats:LinearModel:RankDefDesignMat');