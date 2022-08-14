function out=york_fixed_intercept(data,alpha,afx)

%% A "york fit a with fixed intercept" script modified after an unpublished R script written by Pieter Vermeesch (personal communication)
% This script is similiar to york.m but with a fixed intercept, inspired by
% a R script written by Pieter Vermeesch (personal communication). It's a
% inelegant but quick solution to determine the line.

%% Arguments:
% data must have fields (similar to the scripts in Balco and Rovey, 2008):
%   data.x: x-values (any unit; 1xn vector)
%   data.dx: 1 sigma absolute error of x-values (ditto; 1xn vector)
%   data.y: y-values (ditto; 1xn vector)
%   data.dy: 1 sigma absolute error of y-values (ditto; 1xn vector)
%   data.covxy: correlation coefficients (unitless; 1xn vector; optionally)
% alpha: cutoff value for confidence intervals (unitless; scalar)
% afx: value of the fixed intercept (same as data.x; scalar)

%% Output:
% out have fields as:
%   out.alpha: value of the eponymous input argument (unitless; scalar)
%   out.a: intercept of the straight line fit (same as data.x; scalar)
%   out.sa: standard error of the intercept (ditto; scalar)
%   out.b: slope of the fit (ditto; scalar)
%   out.sb: standard error of the slope (ditto; scalar)
%   out.covab: covariance of the slope and intercept (unitless; scalar)
%   out.df: degrees of freedom of the linear fit (n-2) (unitless; scalar)
%   out.mswd: mean square of the residuals statistic (unitless; scalar)
%   out.pvalue: p-value of a Chi-square value with df degrees of freedom
%   (unitless; scalar) 

    % default values for data.covxy
    n=size(data.x,2); % the number of data
    if (length(fieldnames(data))<5)
        tmp=zeros(1,n);
        data.covxy=tmp;
    end
    a=afx;

    % initial value of the slope
    % The rank of X decreases to 1, because the intercept is constrained.
    X=data.x';
    Y=data.y';
    b=(X'*X)\X'*(Y-afx);

    wX=1./data.dx./data.dx;
    wY=1./data.dy./data.dy;
    A=sqrt(wX.*wY);
    % function squares
    S=@(b)sum(wX.*wY./(wX+b.^2.*wY-2.*b.*data.covxy.*A).*(data.y-b.*data.x-a).^2)/2;
    % optimize the function for least squares and find Hessian matrix
    b0=fminbnd(S,min(b/2,2*b),max(b/2,2*b));
    options = optimoptions('fminunc','Display','off');
    [b,Sval,~,~,~,H]=fminunc(S,b0,options);
    sb=sqrt(inv(H));

    W=wX.*wY./(wX+b.*b.*wY-2.*b.*data.covxy.*A);
    xbar=sum(W.*data.x,'omitnan')/sum(W,'omitnan');
    % output
    out.alpha=alpha;
    out.a=a;
    out.sa=0;
    out.b=b;
    out.sb=sb;
    out.covab=-xbar*sb^2;
    % compute MSWD
    out.df=n-2;
    if (out.df>0)
        out.mswd=Sval*2/out.df;
        out.pvalue=1-chi2cdf(Sval*2,out.df);
    end
    if (out.df<=0)
        out.mswd=1;
        out.pvalue=1;
    end
end
