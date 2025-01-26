function [par, fval] = boundfmin(fun,init,ll,ul,opt)
%BOUNDFMIN Custom function for fitting functions with limited values of
%parameters.

if nargin < 5
    opt = optimset('display', 'off', 'TolX', 1e-5, 'TolFun', 1e-16,'MaxFunEvals',200);
end
% Check there are no mistakes in input limits.
if any(ll > init)
    ll
    init
    error("Requirement: ll < init")
elseif any(init > ul)
    init
    ul
    error("Requirement: ul > init")
elseif any(ll >= ul)
    ll
    ul
    error("Requirement: ll < ul")
end

Delta = ul-ll;

mapfun = @(x) Delta./(1+exp(-x./Delta)) + ll;
invfun = @(x) Delta.*log((x-ll)./(ul-x));

fun2 = @(par) fun(mapfun(par)); % Map inputs from R to [ll,ul]
init = invfun(init); % Map initial values to R
if any(isinf(init))
    % Change infinities to actual values:
    % (infinities will result in errors with fminsearch.)
    init(init==inf) = 1e20;
    init(init==-inf) = realmin;
end

if nargout == 2
    [par,fval] = fminsearch(fun2,init,opt);
elseif nargout ==1
    par = fminsearch(fun2,init,opt);
end

par = mapfun(par); % Map outputs back to specified interval

end