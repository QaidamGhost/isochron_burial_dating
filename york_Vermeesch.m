function out=york_Vermeesch(data,alpha)

% York fit script modified after york.R in IsoplotR by Vermeesch, 2018
%
% "Details: Given n pairs of (approximately) collinear measurements X_i 
% and Y_i (for 1 ≤q i ≤q n), their uncertainties s[X_i] and s[Y_i], and 
% their covariances cov[X_i,Y_i], the york function finds the best fitting 
% straight line using the least-squares algorithm of York et al. (2004). 
% This algorithm is modified from an earlier method developed by 
% York (1968) to be consistent with the maximum likelihood approach of 
% Titterington and Halliday (1979). It computes the MSWD as a measure of 
% under/overdispersion. Overdispersed datasets (MSWD>1) can be dealt with 
% in the same three ways that are described in the documentation of the 
% isochron function."
%
% References: 
% Titterington, D. Michael, and Alex N. Halliday. "On the fitting of 
%   parallel isochrons and the method of maximum likelihood." Chemical 
%   Geology 26.3-4 (1979): 183-195. (https://www.sciencedirect.com/science/
%   article/abs/pii/0009254179900457
% Vermeesch, Pieter. "IsoplotR: A free and open toolbox for geochronology."
%   Geoscience Frontiers 9.5 (2018): 1479-1493. (https://www.sciencedirect.
%   com/science/article/pii/S1674987118300835; https://rdrr.io/cran/
%   IsoplotR/man/york.html
% York, Derek. "Least squares fitting of a straight line with correlated 
%   errors." Earth and planetary science letters 5 (1968): 320-324. (https:
%   //www.sciencedirect.com/science/article/abs/pii/S0012821X68800597
% York, Derek, et al. "Unified equations for the slope, intercept, and 
%   standard errors of the best straight line." American journal of physics
%   72.3 (2004): 367-375. (https://aapt.scitation.org/doi/abs/10.1119/
%   1.1632486
%
% "Arguments:
% data must have fields as:
%   data.x: x-values (nx1 vector)
%   data.dx: absolute error of x-values (nx1 vector)
%   data.y: y-values (nx1 vector)
%   data.dy: absolute error of y-values (nx1 vector)
%   data.covxy: correlation coefficients (nx1 vector; optionally)
% alpha: cutoff value for confidence intervals (scalar; optionally)"
%
% "Output:
% a: intercept of the straight line fit (scalar)
% sa: standard error of the intercept (scalar)
% b: slope of the fit (scalar)
% sb: standard error of the slope (scalar)
% covab: covariance of the slope and intercept (scalar)
% mswd: mean square of the residuals statistic (scalar)
% df: degrees of freedom of the linear fit (n-2) (scalar)
% alpha: value of the eponymous input argument (scalar)
% pvalue: p-value of a Chi-square value with df degrees of freedom (scalar)

    % default values for alpha and data.covxy
    if (nargin<2)
        alpha=0.05;
    end
    n=size(data.x,2); % the number of fields in data
    if (length(fieldnames(data))<5)
        tmp=zeros(1,n);
        data.covxy=tmp;
    end

    % first guess using linear models
    % note that this least square fit method leads to different slopes 
    %   against "stats::lm()" in York.R by Vermeesch
    tmp=data.y/[data.x;ones(1,n)];
    b=tmp(1);

    % main process
    wX=1./data.dx./data.dx;
    wY=1./data.dy./data.dy;
    for i=1:50  % "50 = maximum times for iterations"
        bold=b;
        A=sqrt(wX.*wY);
        W=wX.*wY./(wX+b.*b.*wY-2.*b.*data.covxy.*A);
        Xbar=sum(W.*data.x,'omitnan')/sum(W,'omitnan');
        Ybar=sum(W.*data.y,'omitnan')/sum(W,'omitnan');
        U=data.x-Xbar;
        V=data.y-Ybar;
        B=W.*(U./wY+b.*V./wX-(b.*U+V).*data.covxy/A);
        b=sum(W.*B.*V,'omitnan')/sum(W.*B.*U,'omitnan');
        if ((bold/b-1)^2 < 1e-15)
            break;
        end
    end
    a=Ybar-b*Xbar;
    xbar=sum(W.*data.x,'omitnan')/sum(W,'omitnan');
    u=data.x-xbar;
    sb=sqrt(1/sum(W.*u.*u,'omitnan'));
    sa=sqrt(1/sum(W,'omitnan')+(xbar*sb)^2);

    % output
    out.alpha=alpha;
    out.a=a;
    out.sa=sa;
    out.b=b;
    out.sb=sb;
    out.covab=-xbar*sb^2;
    % compute MSWD
    X2=sum(W.*(data.y-b.*data.x-a).^2);
    out.df=n-2;
    if (out.df>0)
        out.mswd=X2/out.df;
        out.pvalue=1-chi2cdf(X2,out.df);
    end
    if (out.df<=0)
        out.mswd=1;
        out.pvalue=1;
    end
end