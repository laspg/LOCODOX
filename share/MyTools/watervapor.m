function pw=watervapor(T,S)
% function pw=watervapor(T,S)
% calculating pH2O / atm after Weiss and Price 1980
% T in �C
%
% part of optcalc-toolbox
% H. Bittig, IFM-GEOMAR
% 31.03.2011
%
% modification
%   2017 09 28, E. Brion, ALTRAN
%   S and/or T could be a cell array

if iscell(T) & iscell(S)
    pw = cellfun(@(x,y) exp(24.4543-(67.4509*(100./(x+273.15)))-(4.8489*log(((273.15+x)./100)))-0.000544.*y),T,S,'UniformOutput', false);        
elseif iscell(T) & ~iscell(S)
    pw = cellfun(@(x) exp(24.4543-(67.4509*(100./(x+273.15)))-(4.8489*log(((273.15+x)./100)))-0.000544.*S),T,'UniformOutput', false);        
elseif ~iscell(T) & iscell(S)
    pw = cellfun(@(x) exp(24.4543-(67.4509*(100./(T+273.15)))-(4.8489*log(((273.15+T)./100)))-0.000544.*x),S,'UniformOutput', false);        
else
    pw=exp(24.4543-(67.4509*(100./(T+273.15)))-(4.8489*log(((273.15+T)./100)))-0.000544.*S);
end