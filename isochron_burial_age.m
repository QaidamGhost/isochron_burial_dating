function [iso_bur_age,upper_sigma_bur_age,lower_sigma_bur_age] = isochron_burial_age(data,init_Rinh,source_lat,source_elv,measured_lat,measured_elv,shielding_factor,z,rho,option)

%% An iteration process to calculate isochron line for burial dating following Erlanger et al., 2012; Erlanger, 2010; Granger, 2014
%  Note that the script will calculate the upper and lower limits of the  
%  burial age if the intercept of the original isochron line is below zero.
%  A. The upper limit ignores the post-burial production in the
%  samples and calculates the burial age which is similar to simple burial
%  dating.
%  B. On the contrary, the lower limit of burial age is calculated by
%  modelling the post-burial production whose constraint could be
%  different as followings.
%  B1. Gibbon et al., 2009 assumes the depth of their sample has not
%  changed since its burial and erosion has been in absence. This
%  assumption suggests the burial age replacing the exposure duration and
%  calculate the post-burial concentration using eq. 32 in Granger, 2014.
%  Finally, the burial age will be solved iteratively.
%  B2. Unlike this approach, we using a previously known exposure age
%  rather than the burial age to constrain the lower limit. See "Npb_depth.m"
%  for more details of the calculation.

%% Arguments:
% data must have fields as:
%   data.x: measured 10Be concentration (atom/g; 1xn vector)
%   data.dx: 1 sigma absolute error of 10Be (atom/g; 1xn vector)
%   data.y: measured 26Al concentration (atom/g; 1xn vector)
%   data.dy: 1 sigma absolute error of 26Al (atom/g; 1xn vector)
% init_Rinh: initial guess of Rinh (unitless; scalar)
% source_lat: average latitude in the source area (degree; scalar)
% source_elv: average elevation in the source area (m; scalar)
% measured_lat: mearsured latitude of the samples (degree; scalar)
% measured_elv: mearsured elevation of the samples (m; scalar)
% shielding_factor: (unitless; scalar)
% z: depth of samples (cm; scalar)
% rho: density of the overburdens (g/cm^3; scalar)
% option.flag2: "1" for loading e.mat and "0" for using default value zero for
%    erosion rate for a "constant exposure" situation during constraining the
%    lower limit when the first iterated isochron line's intercept is less than
%    zero (unitless; scalar)

%% Output:
% iso_bur_age: isochron burial age (Myr; scalar or 1x3 vector (original
% burial age, upper limit of the burial age, and lower limit of the burial age) 
% if the intercept of the original isochron line is below zero)
% upper_sigma_bur_age: upper 1 sigma absolute error of the burial age
% (ditto)
% lower_sigma_bur_age: lower 1 sigma absolute error of the burial age
% (ditto)

    disp('------------------------------------------------------------------------------------------------------------------');
    disp('Running the program.')
    data_all=data;

    load consts.mat alpha;
    fprintf('Estimated production rate ratio is %.2f.\n',init_Rinh);
    fprintf('The average latitude and elevation in the source area are %.2f degrees and %d meters.\n',source_lat,source_elv);
    if nargin > 4
        fprintf('The latitude and elevation of the samples are %.2f degrees and %d meters.\n',measured_lat,measured_elv);
        fprintf('The shielding factor is %.5f.\n',shielding_factor);
        fprintf('The thickness and density of the overburdens upon the samples are %d cm and %.2f g/cm^3.\n',z,rho);
    end
    fprintf('The cutoff value for confidence intervals is %.2f.\n\n',alpha);
    % simple burial dating
    simple_burial_age(data,source_lat,source_elv,init_Rinh);

    
    % mean life for 10Be, 26Al, and bur
    % 10Be: Chmeleff et al., 2010; Korschinek et al., 2010
    % 26Al: Samworth et al., 1972
    % tau_bur: Granger, 2014, eq. 17
    load consts.mat tau_10 tau_bur sigma_tau_bur;
    load consts.mat limit simulation_times;

    % production rate from spallation for 10Be and 26Al on the surface of the source area
    P100=production_rate(source_lat,source_elv,1,0,2.9,10);
    P260=production_rate(source_lat,source_elv,1,0,2.9,26);
    Rs=P260/P100;
    fprintf('\n');
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
    while true
        out = york(data);
        b=out.b;
        a=out.a;
        if b<0
            disp('Error! Negative slope due to bad data!');
            disp('Try to manually remove outliers or add new 26Al-10Be data that are distinguishable from the former data.');
            iso_bur_age=NaN;
            upper_sigma_bur_age=NaN;
            lower_sigma_bur_age=NaN;
            return;
        end
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
    
    fprintf('MSWD is %.3g.\n',out.mswd);

    % weighted mean for Rinh
    wX=1./data.dx./data.dx;
    wY=1./data.dy./data.dy;
    W=wX.*wY./(wX+b.*b.*wY);
    Rinh_bar=sum(W.*Rinh,'omitnan')/sum(W,'omitnan');

    % Monte-Carlo simulation
    iso_bur_age(1)=-tau_bur*log(b/Rinh_bar);
    cache=zeros(1,simulation_times);
    for i=1:simulation_times
        while true
            rand_b=normrnd(b,sigma_b);
            if rand_b>0
                break;
            end
        end
        rand_tau_bur=normrnd(tau_bur,sigma_tau_bur);
        rand_bur_age=-rand_tau_bur*log(rand_b/Rinh_bar);
        cache(i)=rand_bur_age;
    end
    [~,upper_sigma_bur_age(1),lower_sigma_bur_age(1)]=KDE(cache,iso_bur_age(1));
    %[iso_bur_age(1),upper_sigma_bur_age(1),lower_sigma_bur_age(1)]=KDE(cache);

    plot_isochron(a,sigma_a,b,sigma_b,data,data_backup,removed_data,Rs,option);
    fprintf('Burial age is %f Myr, upper 1 sigma error is %+f Myr, and lower 1 sigma error is %f Myr.\n',iso_bur_age(1),upper_sigma_bur_age(1),lower_sigma_bur_age(1));

    if out.a<0    % intercept less than zero -> use upper and lower limits constraints
        fprintf('\n');
        disp('WARNING! NEGATIVE INTERCEPT DETECTED!');
        if nargin==4    
            disp('PLEASE input measured_lat, measured_elv, shielding_factor, z, rho (, and option) and then run the script again for the upper and lower limits!');
            iso_bur_age=NaN;
            upper_sigma_bur_age=NaN;
            lower_sigma_bur_age=NaN;
            return;
        elseif nargin == 9  % no option for the upper limit
            option.flag2=0;
        end
        disp('The isochron burial age could be dramatically underestimated by the previous value.');
        disp('------------------------------------------------------------------------------------------------------------------');
        disp('Here calculate the upper and lower limits of the burial age.');
        %% upper limit of the burial age: intercept fixed at zero
        disp('Calculating the upper limit of the burial age:');
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
        while true
            out = york_fixed_intercept(data,0);
            b=out.b;
            sigma_b=out.sb;
            if b<0
                disp('Error! Negative slope due to bad data!');
                disp('Try to manually remove outliers or add new 26Al-10Be data that are distinguishable from the former data.');
                iso_bur_age=NaN;
                upper_sigma_bur_age=NaN;
                lower_sigma_bur_age=NaN;
                return;
            end
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
        
        fprintf('MSWD is %.3g.\n',out.mswd);

        % weighted mean for Rinh
        wX=1./data.dx./data.dx;
        wY=1./data.dy./data.dy;
        W=wX.*wY./(wX+b.*b.*wY);
        Rinh_bar=sum(W.*Rinh,'omitnan')/sum(W,'omitnan');

        % Monte-Carlo simulation
        iso_bur_age(2)=-tau_bur*log(b/mean(Rinh_bar));
        cache=zeros(1,simulation_times);
        for i=1:simulation_times
            while true
                rand_b=normrnd(b,sigma_b);
                if rand_b>0
                    break;
                end
            end
            rand_tau_bur=normrnd(tau_bur,sigma_tau_bur);
            rand_bur_age=-rand_tau_bur*log(rand_b/Rinh_bar);
            cache(i)=rand_bur_age;
        end
        [~,upper_sigma_bur_age(2),lower_sigma_bur_age(2)]=KDE(cache,iso_bur_age(2));
        %[iso_bur_age(2),upper_sigma_bur_age(2),lower_sigma_bur_age(2)]=KDE(cache);
        
        plot_isochron(a,sigma_a,b,sigma_b,data,data_backup,removed_data,Rs,option);
        fprintf('Upper limit of the burial age is %f Myr, upper 1 sigma error is %+f Myr, and lower 1 sigma error is %f Myr.\n\n',iso_bur_age(2),upper_sigma_bur_age(2),lower_sigma_bur_age(2));
     
        %% lower limit of the burial age: insert the post-burial concentration as a datum
        disp('Calculating the lower limit of the burial age:');
        data=data_all;
        % the post-burial concentration
        [N10pb,sigma_N10pb,N26pb,sigma_N26pb] = Npb_depth(measured_lat,measured_elv,shielding_factor,z,rho,option);
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
        while true
            out = york(data);
            b=out.b;
            a=out.a;
            sigma_b=out.sb;
            sigma_a=out.sa;
            if b<0
                disp('Error! Negative slope due to bad data!');
                disp('Try to manually remove outliers or add new 26Al-10Be data that are distinguishable from the former data.');
                iso_bur_age=NaN;
                upper_sigma_bur_age=NaN;
                lower_sigma_bur_age=NaN;
                return;
            end
            % recover previous linearized data from data_backup
            data=data_backup;
            % calculate the burial age using slope and Rinh (Erlanger,
            % 2010, eq. 12
            bur_age=-tau_bur*log(b./Rinh);
            % calculate the post burial production of 10Be
            % C10 is determined directly by post-burial 10Be concentration
            C10=data.x(n);
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

        fprintf('MSWD is %.3g.\n',out.mswd);

        % weighted mean for Rinh
        wX=1./data.dx./data.dx;
        wY=1./data.dy./data.dy;
        W=wX.*wY./(wX+b.*b.*wY);
        W(n)=[]; % Rinh of post-burial data point is precluded
        Rinh_bar=sum(W.*Rinh,'omitnan')/sum(W,'omitnan');

        % Monte-Carlo simulation
        iso_bur_age(3)=-tau_bur*log(b/Rinh_bar);
        cache=zeros(1,simulation_times);
        for i=1:simulation_times
            while true
                rand_b=normrnd(b,sigma_b);
                if rand_b>0
                    break;
                end
            end
            rand_tau_bur=normrnd(tau_bur,sigma_tau_bur);
            rand_bur_age=-rand_tau_bur*log(rand_b/Rinh_bar);
            cache(i)=rand_bur_age;
        end
        [~,upper_sigma_bur_age(3),lower_sigma_bur_age(3)]=KDE(cache,iso_bur_age(3));
        %[iso_bur_age(3),upper_sigma_bur_age(3),lower_sigma_bur_age(3)]=KDE(cache);

        [P10n_z,~,P10ms_z,~,P10mf_z,~]=production_rate(measured_lat,measured_elv,shielding_factor,z,rho,10);
        [P26n_z,~,P26ms_z,~,P26mf_z,~]=production_rate(measured_lat,measured_elv,shielding_factor,z,rho,26);
        option.Rd=(P26n_z+P26ms_z+P26mf_z)/(P10n_z+P10ms_z+P10mf_z);
        plot_isochron(a,sigma_a,b,sigma_b,data,data_backup,removed_data,Rs,option);
        fprintf('Lower limit of the burial age is %f Myr, upper 1 sigma error is %+f Myr, and lower 1 sigma error is %f Myr.\n',iso_bur_age(3),upper_sigma_bur_age(3),lower_sigma_bur_age(3));
    end
    disp('------------------------------------------------------------------------------------------------------------------');
end
