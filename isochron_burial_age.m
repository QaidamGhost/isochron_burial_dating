function [iso_bur_age,upper_sigma_bur_age,lower_sigma_bur_age] = isochron_burial_age(data,init_Rinh,limit,source_lat,source_elv,measured_lat,measured_elv,z,rho,alpha)

%% An iteration process to calculate isochron line for burial dating following Erlanger et al., 2012; Erlanger, 2010; Granger, 2014
%  Note that the script will calculate minimum and maximum burial age if
%  the intercept of the original isochron line is below zero. 
%  A. The minimum estimation ignores the post-burial production in the
%  samples and calculates the burial age which is similar to simple burial
%  dating.
%  B. On the contrary, the maximum estimation calculate the burial age by
%  maximizing the post-burial production whose constraint could be
%  different as followings.
%  B1. Gibbon et al., 2009 assumes the depth of their sample has not
%  changed since its burial and erosion has been in absence. This
%  assumption suggests the burial age replacing the exposure duration and
%  calculate the post-burial concentration using eq. 32 in Granger, 2014.
%  Finally, the burial age will be solved iteratively.
%  B2. Unlike this approach, we using a previously known exposure age
%  rather than the burial age to estimate the maximum. See "Npb_depth.m"
%  for more details of the calculation.

%% Arguments:
% data must have fields as:
%   data.x: measured 10Be concentration (atom/g; 1xn vector)
%   data.dx: 1 sigma absolute error of 10Be (atom/g; 1xn vector)
%   data.y: measured 26Al concentration (atom/g; 1xn vector)
%   data.dy: 1 sigma absolute error of 26Al (atom/g; 1xn vector)
% init_Rinh: initial guess of Rinh (unitless; scalar)
% limit: iteration stops if the variation of slope reaches the given limit
% (unitless; scalar)
% source_lat: average latitude in the source area (degree; scalar)
% source_elv: average elevation in the source area (m; scalar)
% measured_lat: mearsured latitude of the samples (degree; scalar)
% measured_elv: mearsured elevation of the samples (m; scalar)
% z: depth of samples (cm; scalar)
% rho: density of the overburdens (g/cm^3; scalar)
% alpha: cutoff value for confidence intervals (unitless; scalar)

%% Output:
% iso_bur_age: isochron burial age (Myr; scalar or 1x3 vector (original
% burial age, minimum burial age, and maximum burial age) if the intercept
% of the original isochron line is below zero)
% upper_sigma_bur_age: upper 1 sigma absolute error of the burial age
% (ditto)
% lower_sigma_bur_age: lower 1 sigma absolute error of the burial age
% (ditto)

    disp('------------------------------------------------------------------------------------------------------------------');
    disp('Running the program.')
    data_all=data;

    % mean life for 10Be, 26Al, and bur
    tau_10=2.001;   % Chmeleff et al., 2010; Korshchinek et al., 2010
    sigma_tau_10=0.017;
    tau_26=1.034;   % Samworth et al., 1972
    sigma_tau_26=0.024;
    tau_bur=1/(1/tau_26-1/tau_10);  % Granger, 2014, eq. 17
    sigma_tau_bur=sqrt((tau_10^2/(tau_10-tau_26)^2)^2*sigma_tau_26^2+(tau_26^2/(tau_10-tau_26)^2)^2*sigma_tau_10^2);

    % production rate for 10Be and 26Al on the surface of the source area
    [Pn,~,Pms,~,Pmf,~]=production_rate(source_lat,source_elv,0,rho,10);
    P100=Pn+Pms+Pmf;
    [Pn,~,Pms,~,Pmf,~]=production_rate(source_lat,source_elv,0,2.65,26);
    P260=Pn+Pms+Pmf;
    disp('Calculating the isochron burial age:');
    % preclude reworked samples
    option.flag=0;
    [data,removed_data,~] = remove_outliers(data,alpha,option);

    %% isochron burial dating
    % initialization
    data_backup=data;
    count=0;    % iteration times
    old_b=0;    % slope determined by previous iteration
    n=size(data.x,2);   % number of samples
    if n==1
        disp('Error! Only one sample!');
        iso_bur_age=NaN;
        upper_sigma_bur_age=NaN;
        lower_sigma_bur_age=NaN;
        return;
    end
    if n==2
        disp('Only two samples! Isochron-line is fully determined!');
    end
    Rinh=zeros(1,n);
    for i=1:n
        Rinh(1,i)=init_Rinh;    % Rinh for each sample
    end

    % iteration start
    while(1)
        out = york(data,alpha);
        b=out.b;
        a=out.a;
        sigma_a=out.sa;
        sigma_b=out.sb;
        % recover previous linearized data from data_backup
        data=data_backup;
        % calculate the burial age using slope and Rinh (Erlanger, 2010,
        % eq. 12
        bur_age=-tau_bur*log(b./Rinh);
        % calculate the post burial production of 10Be
        C10=a/(P260/P100-b);
        % calculate the inherited concentration of 10Be for each sample
        N10inh=(data.x-C10).*exp(bur_age./tau_10);
        % calculate the inheritance ratio of 26Al/10Be for each sample
        Rinh=P260./(P100+N10inh./tau_bur/1E6);
        % calculate the linearization factor for each sample (Erlanger,
        % 2012
        LF=Rinh./(P260/P100);
        % linearize 10Be for the next regression
        data.x=data.x.*LF;
        data.dx=data.dx.*LF;
        % proceed unconditionally during the first iteration
        if count>0
            % iteration breaks after 100 attempts
            if (abs(old_b-b)<limit) || (count>=99)
               break;
            end
        end
        old_b=b;
        count=count+1;
    end
    
    % Monte-Carlo simulation
    cache=zeros(1,1E5);
    for i=1:1E5
        rand_b=normrnd(b,sigma_b);
        rand_tau_bur=normrnd(tau_bur,sigma_tau_bur);
        rand_bur_age=-rand_tau_bur*log(rand_b/mean(Rinh));
        cache(i)=rand_bur_age;
    end
    [iso_bur_age(1),upper_sigma_bur_age(1),lower_sigma_bur_age(1)]=KDE(cache);

    plot_isochron(a,sigma_a,b,sigma_b,data,data_backup,removed_data,init_Rinh,option);
    fprintf('Burial age is %f Myr, upper 1 sigma error is %f Myr, and lower 1 sigma error is %f Myr.\n',iso_bur_age(1),upper_sigma_bur_age(1),lower_sigma_bur_age(1));

    if out.a<0    % intercept less than zero -> use max and min estimation
        fprintf('\n');
        disp('WARNING! NEGATIVE INTERCEPT DETECTED!');
        disp('The isochron burial age could be dramatically underestimated by the previous value.');
        disp('------------------------------------------------------------------------------------------------------------------');
        disp('Here calculate the minimum and maximum constraints on the burial age.');
        %% min estimation: intercept fixed at zero
        disp('Calculating the minimum burial age:');
        data=data_all;
        option.flag=1;
        [data,removed_data,~] = remove_outliers(data,alpha,option);
        % initialization
        data_backup=data;
        count=0;    % iteration times
        old_b=0;    % slope determined by previous iteration
        n=size(data.x,2);   % number of samples
        if n==1
            disp('Error! Only one sample!');
            iso_bur_age=NaN;
            upper_sigma_bur_age=NaN;
            lower_sigma_bur_age=NaN;
            return;
        end
        if n==2
            disp('Only two samples! Isochron-line is fully determined!');
        end
        Rinh=zeros(1,n);
        for i=1:n
            Rinh(1,i)=init_Rinh;    % Rinh for each sample
        end
        a=0;
        sigma_a=0;

        % iteration start
        while(1)
            out = york_fixed_intercept(data,alpha,0);
            b=out.b;
            sigma_b=out.sb;
            % recover previous linearized data from data_backup
            data=data_backup;
            % calculate the burial age using slope and Rinh (Erlanger,
            % 2010, eq. 12 
            bur_age=-tau_bur*log(b./Rinh);
            % calculate the inherited concentration of 10Be for each sample
            N10inh=data.x.*exp(bur_age./tau_10);
            % calculate the inheritance ratio of 26Al/10Be for each sample
            Rinh=P260./(P100+N10inh./tau_bur/1E6);
            % calculate the linearization factor for each sample (Erlanger,
            % 2012
            LF=Rinh./(P260/P100);
            % linearize 10Be for the next regression
            data.x=data.x.*LF;
            data.dx=data.dx.*LF;
            % proceed unconditionally during the first iteration
            if count>0
                % iteration breaks after 100 attempts
                if (abs(old_b-b)<limit) || (count>=99)
                   break;
                end
            end
            old_b=b;
            count=count+1;
        end
        
        % Monte-Carlo simulation
        % If your final minimum burial age is unstable after several
        % repeats, try to set a larger value for "num" so that the script
        % will choose a median value from the simulations. However, this
        % will take longer time to get the final result.
        num=1;
        for j=1:(2*num+1)
	        cache=zeros(1,1E5);
            for i=1:1E5
                rand_b=normrnd(b,sigma_b);
                rand_tau_bur=normrnd(tau_bur,sigma_tau_bur);
                rand_bur_age=-rand_tau_bur*log(rand_b/mean(Rinh));
                cache(i)=rand_bur_age;
            end
            [most_prob,upper_sigma,lower_sigma]=KDE(cache);
            rand.bur_age(j)=most_prob;
            rand.upper_sigma_bur_age(j)=upper_sigma;
            rand.lower_sigma_bur_age(j)=lower_sigma;
        end
        m=median(rand.bur_age);
        for j=1:(2*num+1)
            if rand.bur_age(j)==m
                break;
            end
        end
        iso_bur_age(2)=m;
        upper_sigma_bur_age(2)=rand.upper_sigma_bur_age(j);
        lower_sigma_bur_age(2)=rand.lower_sigma_bur_age(j);
        
        %
        plot_isochron(a,sigma_a,b,sigma_b,data,data_backup,removed_data,init_Rinh,option);
        fprintf('Minimum burial age is %f Myr, upper 1 sigma error is %f Myr, and lower 1 sigma error is %f Myr.\n\n',iso_bur_age(2),upper_sigma_bur_age(2),lower_sigma_bur_age(2));
        
        %% max estimation: insert the post-burial concentration as a datum
        disp('Calculating the maximum burial age:');
        data=data_all;
        option.flag2=0;
        % the post-burial concentration
        [N10pb,sigma_N10pb,N26pb,sigma_N26pb] = Npb_depth(measured_lat,measured_elv,z,rho,option);
        option.flag=2;
        option.Npb.x=N10pb;
        option.Npb.dx=sigma_N10pb;
        option.Npb.y=N26pb;
        option.Npb.dy=sigma_N26pb;
        [data,removed_data,~] = remove_outliers(data,alpha,option);
        % combine the post-burial concentration with the measured data
        % the post-burial concentration is the last elements in each field
        % in data
        data.x=[data.x,N10pb];
        data.dx=[data.dx,sigma_N10pb];
        data.y=[data.y,N26pb];
        data.dy=[data.dy,sigma_N26pb];
        % initialization
        data_backup=data;
        count=0;    % iteration times
        old_b=0;    % slope determined by previous iteration
        n=size(data.x,2);   % number of samples %!
        if n==1
            disp('Error! Only one sample!');
            iso_bur_age=NaN;
            upper_sigma_bur_age=NaN;
            lower_sigma_bur_age=NaN;
            return;
        end
        if n==2
            disp('Only two samples! Isochron-line is fully determined!');
        end
        Rinh=zeros(1,n);
        for i=1:n
            Rinh(1,i)=init_Rinh;    % Rinh for each sample
        end

        % iteration start
        while(1)
            out = york(data,alpha);
            b=out.b;
            a=out.a;
            sigma_b=out.sb;
            sigma_a=out.sa;
            % recover previous linearized data from data_backup
            data=data_backup;
            % calculate the burial age using slope and Rinh (Erlanger,
            % 2010, eq. 12
            bur_age=-tau_bur*log(b./Rinh);
            % calculate the post burial production of 10Be
            C10=a/(P260/P100-b);
            % calculate the inherited concentration of 10Be for each sample
            N10inh=(data.x-C10).*exp(bur_age./tau_10);
            % calculate the inheritance ratio of 26Al/10Be for each sample
            Rinh=P260./(P100+N10inh./tau_bur/1E6);
            % calculate the linearization factor for each sample (Erlanger,
            % 2012
            LF=Rinh./(P260/P100);
            % The post_burial concentration should not be linearized,
            % because its inheritance component is zero.
            LF(n)=1; 
            % linearize 10Be for the next regression
            data.x=data.x.*LF;
            data.dx=data.dx.*LF;
            % proceed unconditionally during the first iteration
            if count>0
                % iteration breaks after 100 attempts
                if (abs(old_b-b)<limit) || (count>=99)
                   break;
                end
            end
            old_b=b;
            count=count+1;
        end
        Rinh(n)=[]; % remove the post-burial's value from Rinh
        
        % Monte-Carlo simulation
        % If your final maximum burial age is unstable after several
        % repeats, try to set a larger value for "num" so that the script
        % will choose a median value from the simulations. However, this
        % will take longer time to get the final result.
        num=5;
        for j=1:(2*num+1)
            cache=zeros(1,1E5);
            for i=1:1E5
                rand_b=normrnd(b,sigma_b);
                rand_tau_bur=normrnd(tau_bur,sigma_tau_bur);
                rand_bur_age=-rand_tau_bur*log(rand_b/mean(Rinh));
                cache(i)=rand_bur_age;
            end
            [most_prob,upper_sigma,lower_sigma]=KDE(cache);
            rand.bur_age(j)=most_prob;
            rand.upper_sigma_bur_age(j)=upper_sigma;
            rand.lower_sigma_bur_age(j)=lower_sigma;
        end
        m=median(rand.bur_age);
        for j=1:(2*num+1)
            if rand.bur_age(j)==m
                break;
            end
        end
        iso_bur_age(3)=m;
        upper_sigma_bur_age(3)=rand.upper_sigma_bur_age(j);
        lower_sigma_bur_age(3)=rand.lower_sigma_bur_age(j);

        plot_isochron(a,sigma_a,b,sigma_b,data,data_backup,removed_data,init_Rinh,option);
        fprintf('Maximum burial age is %f Myr, upper 1 sigma error is %f Myr, and lower 1 sigma error is %f Myr.\n',iso_bur_age(3),upper_sigma_bur_age(3),lower_sigma_bur_age(3));
        disp('------------------------------------------------------------------------------------------------------------------');
    end
end
