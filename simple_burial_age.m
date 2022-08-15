function [simple_bur_age,upper_sigma_bur_age,lower_sigma_bur_age] = simple_burial_age(data,source_lat,source_elv,limit,init_Rinh)

%% Calculate the simple burial age of each sample

%% Arguments:
% data must have fields as:
%   data.x: measured 10Be concentration (atom/g; 1xn vector)
%   data.dx: 1 sigma absolute error of 10Be (atom/g; 1xn vector)
%   data.y: measured 26Al concentration (atom/g; 1xn vector)
%   data.dy: 1 sigma absolute error of 26Al (atom/g; 1xn vector)
% source_lat: average latitude in the source area (degree; scalar)
% source_elv: average elevation in the source area (m; scalar)

%% Output:
% simple_bur_age: simple burial age (Myr; 1xn vector)
% upper_sigma_bur_age: upper 1 simga absolute error of the burial age (Myr;
% 1xn vector) 
% lower_sigma_bur_age: lower 1 simga absolute error of the burial age (Myr;
% 1xn vector) 

    disp('Calculating the simple burial age:');
    tau_10=2.001;   % Chmeleff et al., 2010; Korshchinek et al., 2010
    sigma_tau_10=0.017;
    tau_26=1.034;   % Samworth et al., 1972
    sigma_tau_26=0.024;
    tau_bur=1/(1/tau_26-1/tau_10);  % Granger, 2014, eq. 17
    sigma_tau_bur=sqrt((tau_10^2/(tau_10-tau_26)^2)^2*sigma_tau_26^2+(tau_26^2/(tau_10-tau_26)^2)^2*sigma_tau_10^2);
    
    % defines of the variables
    N10=data.x;
    sigma_N10=data.dx;
    N26=data.y;
    sigma_N26=data.dy;
    R=N26./N10;
    P100=production_rate(source_lat,source_elv,1,0,2.9,10);
    P260=production_rate(source_lat,source_elv,1,0,2.9,26);

    n=size(N10,2);   % number of data
    limit=ones(1,n)*limit;
    count=0;
    t=zeros(1,n);
    old_t=zeros(1,n);
    Rinh=ones(1,n)*init_Rinh;
    for i=1:n
        while true
            t(i) = -tau_bur*log(R(i)/Rinh(i));
            if count>0
                if abs(old_t(i)-t(i))<limit
                    break;
                end
                if count>=99
                    disp('Error! No iterative solution found!');
                    simple_bur_age=NaN;
                    upper_sigma_bur_age=NaN;
                    lower_sigma_bur_age=NaN;
                    return;
                end
            end
            Rinh(i)=P260/(P100+N10(i)*exp(t(i)/tau_10)/tau_bur/1E6);
            count=count+1;
            old_t(i)=t(i);
        end
    end

    simple_bur_age=t;
    cache=zeros(1E5,n); % malloc
    for i=1:1E5
        rand_N10=normrnd(N10,sigma_N10);
        rand_N26=normrnd(N26,sigma_N26);
        rand_tau_bur=normrnd(tau_bur,sigma_tau_bur);
        rand_burial_age=-rand_tau_bur.*log(rand_N26./rand_N10./Rinh);
        cache(i,:)=rand_burial_age;
    end
    upper_sigma_bur_age=zeros(1,n);
    lower_sigma_bur_age=zeros(1,n);
    for i=1:n
        single=cache(:,i);
        [~,upper_sigma_bur_age(i),lower_sigma_bur_age(i)]=KDE(single);
    end
    fprintf('Index  |  measured 10Be (at/g)  |  1 sigma error (at/g)  |  measured 26Al (at/g)  |  1 sigma error (at/g)  |  Simple burial age (Myr)  |  Upper 1 sigma error (Myr)  |  Lower 1 sigma error (Myr)\n');
    for i=1:n
        fprintf('%5d  |  %20d  |  %20d  |  %20d  |  %20d  |  %23.3g  |  %+25.3g  |  %25.3g\n',i,N10(i),sigma_N10(i),N26(i),sigma_N26(i),simple_bur_age(i),upper_sigma_bur_age(i),lower_sigma_bur_age(i));
    end
end
