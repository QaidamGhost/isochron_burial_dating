function [bur_age,sigma_bur_age] = isochron_Erlanger(data,init_Rinh, ...
    limit,P100,P260)

% An iteration process to calculate isochron line for burial dating 
% following Erlanger et al., 2012; Erlanger, 2010; Granger, 2014
% 
% References:
% Erlanger, Erica D., Darryl E. Granger, and Ryan J. Gibbon. "Rock uplift 
%   rates in South Africa from isochron burial dating of fluvial and marine 
%   terraces." Geology 40.11 (2012): 1019-1022. (https://pubs.
%   geoscienceworld.org/gsa/geology/article/40/11/1019/130729/Rock-uplift
%   -rates-in-South-Africa-from-isochron
% Erlanger, Erica D. "Rock uplift, erosion, and tectonic uplift of South
%   Africa determined with cosmogenic aluminum-26 and beryllium-10." Ph. D.
%   Thesis (2010). (https://www.proquest.com/docview/861934022/
%   89546D55724C480APQ/1
% Granger, D. E. "Cosmogenic nuclide burial dating in archaeology and 
%   paleoanthropology." (2014): 81-97. (https://www.sciencedirect.com/
%   science/article/pii/B9780080959757012080
%
% Arguments:
% data must have fields as:
%   data.x: measured 10Be concentration (atom/g; nx1 vector)
%   data.dx: absolute error of 10Be (atom/g; nx1 vector)
%   data.y: measured 26Al concentration (atom/g; nx1 vector)
%   data.dy: absolute error of 26Al (atom/g; nx1 vector)
% init_Rinh: initial guess of Rinh (unitless; scalar)
% P100: 10Be production rate at surface in source area (atom/g/yr; scalar)
% P260: 26Al production rate at surface in source area (atom/g/yr; scalar)
% limit: iteration stops if the variation of slope reaches the given limit 
%   (unitless; scalar)
% 
% Output:
% bur_age: isochron burial age (Myr; scalar)
% sigma_bur_age: error of burial age (Myr; scalar)
%
% Example:
% Yang:
% data.x=[685595,1341794,867201,1111673,1219183]
% data.dx=[15457,32532,20433,21054,35624]
% data.y=[3117757,4409884,4053256,5554372,5733781]
% data.dy=[86074,137590,130328,158326,179906]
% init_Rinh=6.61
% P260=425.9379
% P100=63.9434
% limit=1E-12

    % mean life for 10Be, 26Al, and bur
    tau_10=2.001;   % Chmeleff et al., 2010; Korshchinek et al., 2010
    % sigma_tau_10=0.017;
    tau_26=1.017;   % Norris et al., 1983
    % sigma_tau_26=0.035;
    tau_bur=1/(1/tau_26-1/tau_10);  % Granger, 2014, eqn17
    
    % preclude reworked samples
    [data,removed_data,~] = kill_outliers_Odom(data,.05);
    % initialization
    data_backup=data;
    count=0;    % iteration times
    old_b=0;    % slope determined by previous iteration
    n=size(data.x,2);   % number of samples
    if n==1
        disp('Error! Only one sample!');
        bur_age=NaN;
        sigma_bur_age=NaN;
        return;
    end
    if n==2
        disp('Only two samples! Isochron-line is fully determined!')
    end
    Rinh=zeros(1,n);
    for i=1:n
        Rinh(1,i)=init_Rinh;    % Rinh for each sample
    end

    % iteration start
    while(1)
        out = york_Vermeesch(data,0.05);
        b=out.b;
        a=out.a;
        % recover previous linearized data from data_backup
        data=data_backup;
        % calculate the burial age using slope and Rinh (Erlanger, 2010,...
        %   eqn12
        bur_age=-tau_bur*log(b./Rinh);
        % calculate the post burial production of 10Be
        C10=a/(P260/P100-b);
        % calculate the inherited concentration of 10Be for each sample
        N10inh=(data.x-C10).*exp(bur_age./tau_10);
        % calculate the inheritance ratio of 26Al/10Be for each sample
        Rinh=P260./(P100+N10inh./tau_bur/1E6);
        % calculate the linearization factor for each sample (Erlanger,...
        %   2012
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

    % final caculation of the burial age
    bur_age=-tau_bur*log(b/mean(Rinh));
    % re-caculate the error of the converged slope
    sigma_b=sigma_slope(b,data);
    % monte-carlo simulation for 100000 times
    cache=zeros(1,1E5);
    for i=1:1E5
        rand_b=normrnd(b,sigma_b);
        rand_bur_age=-tau_bur*log(rand_b/mean(Rinh));
        cache(i)=rand_bur_age;
    end
    % error of the burial age
    sigma_bur_age=std(cache);
    % plot
    plot_isochron(a,b,sigma_b,data,data_backup,removed_data,init_Rinh);
    fprintf('Burial age is %f Myr, 1 sigma error is %f Myr.\n',bur_age,sigma_bur_age);
end
