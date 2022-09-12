function out=york(data)

%% A "york fit" script modified after york.R in IsoplotR by Vermeesch, 2018
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

%% Arguments:
% data must have fields (similar to the scripts in Balco and Rovey, 2008):
%   data.x: x-values (any unit; 1xn vector)
%   data.dx: 1 sigma absolute error of x-values (ditto; 1xn vector)
%   data.y: y-values (ditto; 1xn vector)
%   data.dy: 1 sigma absolute error of y-values (ditto; 1xn vector)
%   data.covxy: correlation coefficients (unitless; 1xn vector; optionally)

%% Output:
% out have fields as:
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

    % first guess using linear models
    % Note that this least square fit method leads to different initial
    % value of slopes against "stats::lm()" in york.R by Vermeesch. Yet,
    % the value of convergent slope will be the same using both scripts.
    X=[data.x;ones(1,n)]';
    Y=data.y';
    tmp=(X'*X)\X'*Y;
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
        B=W.*(U./wY+b.*V./wX-(b.*U+V).*data.covxy./A);
        b=sum(W.*B.*V,'omitnan')/sum(W.*B.*U,'omitnan');
        if ((bold/b-1)^2 < 1e-15)
            break;
        end
    end
    a=Ybar-b*Xbar;
    xbar=sum(W.*data.x,'omitnan')/sum(W,'omitnan');
    u=data.x-xbar;
    % calculation of 1 sigma error of the slope and intercept
    % The two following equations are equivalent to eq. 17 and 18,
    % substituted by eq. from 19, 20, 21, 22, and 23, in Mahon, 1996.
    sb=sqrt(1/sum(W.*u.*u,'omitnan'));
    sa=sqrt(1/sum(W,'omitnan')+(xbar*sb)^2);
   
    % output
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
