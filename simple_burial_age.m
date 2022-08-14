function [bur_age,upper_sigma_bur_age,lower_sigma_bur_age]=simple_burial_age(data,measured_lat,measured_elv)

%% Calculate the simple burial age of each sample following the equation 22 in Granger, 2014.
% Note that the formula of the burial age do not give an analytic solution
% but an simple estimation. See "Granger, 2014" for more imformation.

%% Arguments:
% data must have fields as:
%   data.x: measured 10Be concentration (atom/g; 1xn vector)
%   data.dx: 1 sigma absolute error of 10Be (atom/g; 1xn vector)
%   data.y: measured 26Al concentration (atom/g; 1xn vector)
%   data.dy: 1 sigma absolute error of 26Al (atom/g; 1xn vector)
% measured_lat: mearsured latitude of the samples (degree; scalar)
% measured_elv: mearsured elevation of the samples (m; scalar)

%% Output:
% bur_age: isochron burial age (Myr; 1xn vector)
% upper_sigma_bur_age: upper 1 simga absolute error of the burial age (Myr;
% 1xn vector) 
% lower_sigma_bur_age: lower 1 simga absolute error of the burial age (Myr;
% 1xn vector) 

    m=size(measured_lat,2);   % number of measured_lat
    % defines of the variables
    N10=data.x;
    sigma_N10=data.dx;
    N26=data.y;
    sigma_N26=data.dy;
    % malloc
    P100=zeros(1,m);
    sigma_P100=P100;
    P260=P100;
    sigma_P260=P100;
    for i=1:m
        [Pn,sigma_Pn,Pms,sigma_Pms,Pmf,sigma_Pmf]=production_rate(measured_lat(i),measured_elv(i),0,2.65,10);
        P100(i)=Pn+Pms+Pmf;
        sigma_P100(i)=sigma_Pn+sigma_Pms+sigma_Pmf;
        [Pn,sigma_Pn,Pms,sigma_Pms,Pmf,sigma_Pmf]=production_rate(measured_lat(i),measured_elv(i),0,2.65,26);
        P260(i)=Pn+Pms+Pmf;
        sigma_P260(i)=sigma_Pn+sigma_Pms+sigma_Pmf;
    end

    % mean life for 10Be and 26Al
    tau_10=2.001;   % Chmeleff et al., 2010; Korshchinek et al., 2010
    sigma_tau_10=0.017;
    tau_26=1.034;   % Samworth et al., 1972
    sigma_tau_26=0.024;

    % Monte-Carlo simulation for the 1 sigma error of the burial age
    n=size(data.x,2);   % number of samples
    cache=zeros(1E5,n); % malloc
    % The errors of the mean life, measured concentrations, production rate
    % at the local surface of the 10Be and 26Al are all taken into the
    % error of the simple burial age.
    for i=1:1E5
        rand_tau_10=normrnd(tau_10,sigma_tau_10);
        rand_tau_26=normrnd(tau_26,sigma_tau_26);
        rand_N10=normrnd(N10,sigma_N10);
        rand_N26=normrnd(N26,sigma_N26);
        rand_P100=normrnd(P100,sigma_P100);
        rand_P260=normrnd(P260,sigma_P260);
        rand_age=(rand_tau_10+1/(1/rand_tau_26-1/rand_tau_10))/2.*log(-0.5./rand_N10.*rand_P100.*(rand_tau_10*1E+06)+sqrt(0.25./(rand_N10.^2).*(rand_P100.^2).*((rand_tau_10*1E+06)^2)+2./rand_N26.*rand_P260.*(rand_tau_26*1E+06)));   % Granger, 2015, eq. 22
        cache(i,:)=rand_age;
    end
    % malloc
    bur_age=zeros(1,n);
    upper_sigma_bur_age=zeros(1,n);
    lower_sigma_bur_age=zeros(1,n);
    for i=1:n
        single=cache(:,i);
        [bur_age(i),upper_sigma_bur_age(i),lower_sigma_bur_age(i)]=KDE(single);
    end
end
