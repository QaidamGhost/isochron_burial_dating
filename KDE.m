function [most_prob,upper_sigma,lower_sigma] = KDE(cache)

%% Find the most probable value and 1 sigma absolute error for the observed simulations of an random variable after establishing the probability density function by a "nonparametric estimation" approach knows as kernel density estimation (KDE).

%% Arguments: 
% cache: the observed simulations of an random variable to be evaluated
% (1×n vector)

%% Output:
% most_prob: the most probable value of the variable (scalar)
% upper_sigma: the upper 1 sigma absolute error of the variable (scalar)
% lower_sigma: the lower 1 sigma absolute error of the variable (scalar)

    step=1E-4;
    min_cache=min(cache);
    max_cache=max(cache);
    bounded_cache=min_cache:step:max_cache;
    % search the most probable value from the simulations
    f=ksdensity(cache,bounded_cache);
    [~,index]=max(f);
    most_prob=bounded_cache(index);
    % fake-normalization with the given step (F to 1)
    nf=f/sum(f*step);
    most_prob_cdp=sum(nf(1:index))*step;
    p=normcdf([-1 0]);
    p1sigma=p(2)-p(1);
    upper_cdp=most_prob_cdp+p1sigma;
    lower_cdp=most_prob_cdp-p1sigma;
    % search for the +-1 sigma value
    n=max(size(f));
    nF=zeros(1,n);
    upper_sigma=0;
    lower_sigma=0;
    % cumulation of the probability by calculus
    for j=1:n
        if j==1
            nF(j)=nF(j)+nf(j)*step;
        elseif j~=1
            nF(j)=nF(j-1)+nf(j)*step;
        end
        if nF(j)>=upper_cdp&&upper_sigma==0
            upper_sigma=min_cache+step*(j-1);
            upper_sigma=upper_sigma-most_prob;
        end
        if nF(j)>=lower_cdp&&lower_sigma==0
            lower_sigma=min_cache+step*(j-1);
            lower_sigma=lower_sigma-most_prob;
        end
    end
end