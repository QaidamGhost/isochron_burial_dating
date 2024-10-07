function [N10pb,sigma_N10pb,N26pb,sigma_N26pb] = Npb_depth(measured_lat,measured_elv,shielding_factor,z,rho,option)

%% Calculate the upper limit of the post-burial concentrations with 1 sigma error at the sampling location and penetration depth if the exposure age (and additional erosion rate) of the surface is constrained.
%
%% Arguments:
% measured_lat: mearsured latitude of the samples (degree; scalar)
% measured_elv: mearsured elevation of the samples (m; scalar)
% shielding_factor: (unitless; scalar)
% z: depth of the samples (cm; scalar)
% rho: density of the overburdens (g/cm^3; scalar)
% option have a field as:
%   option.flag2: "1" for loading e.mat and "0" for using default value
%   zero for erosion rate for a "constant exposure" situation (unitless;
%   scalar)
%
%% Output
% N10pb: the upper limit of the post burial 10Be concentration
% (atom/g;scalar) 
% sigma_N10pb: 1 sigma absolute error of N10pb (atom/g;scalar)
% N26pb: the upper limit of the post burial 26Al concentration
% (atom/g;scalar) 
% sigma_N26pb: 1 sigma absolute error of N26pb (atom/g;scalar)

% The following two .mat files are needed in the calculation. They are the
% results of simulations from a MATLAB program "10Be_profile_simulator" by
% Hidy et al., 2010, G3. Contact Yizhou Yang (yyz606@pku.edu.cn) if you
% have troubles in exporting the simulated results.
% A.The "expo_age.mat" refers to exposure age. It has two variables
% "expo_age" (Kyr; scalar; the most probable value of the simulations of
% exposure age) and "expo_age_est" (Kyr; nx1 vector; the simulations of
% exposure age during each simulation).
% B. The "e.mat" refers to the erosion rate at the surface. It has two
% variables "expo_age" (cm/Kyr; scalar; the most probable value of the
% simulations of erosion rate) and "expo_age_est" (cm/Kyr; nx1 vector; the
% simulations of erosion rate during each simulation).
% If you want to calculate upper limit of the isochron burial age and 
% already know the value and absolate 1 sigma error of erosion rate and 
% exposure age, you can use "custom_mat_files.m" to generate "e.mat" and 
% "expo_age.mat".
% Note that the input exposure age and erosion rate should be Gaussian
% distribution, or you should use Hidy's modified scripts to generate the
% "expo_age.mat" and "e.mat" file.

    if option.flag2==1
        load e.mat e e_est
        e=e/1E+03;
        e_est=e_est/1E+03;
        n2=max(size(e_est));
    elseif option.flag2==0
        e=0;
    end
    load expo_age.mat expo_age expo_age_est
    expo_age_est=expo_age_est./1E3;
    expo_age=expo_age/1E3;
    n=max(size(expo_age_est));

    % mean life for 10Be, 26Al, and bur
    % 10Be: Chmeleff et al., 2010; Korschinek et al., 2010
    % 26Al: Samworth et al., 1972
    % tau_bur: Granger, 2014, eq. 17
    load consts.mat tau_10 sigma_tau_10 tau_26 sigma_tau_26;
    load consts.mat simulation_times;

    [P10n_z,sigma_P10n_z,P10ms_z,sigma_P10ms_z,P10mf_z,sigma_P10mf_z]=production_rate(measured_lat,measured_elv,shielding_factor,z,rho,10);
    [P26n_z,sigma_P26n_z,P26ms_z,sigma_P26ms_z,P26mf_z,sigma_P26mf_z]=production_rate(measured_lat,measured_elv,shielding_factor,z,rho,26);

    %% Attenuation lengths of neutrons, slow muons, and fast muons, respectively. 
    % Unit: g/cm^2. Braucher et al., 2011.
    load consts.mat Ln Lms Lmf;
    
    % post-burial 10Be (eqn6 in Lal, 1991)
    N10n_pb=P10n_z/(1/tau_10/1E+06+rho*e/Ln)*(1-exp(-expo_age/tau_10-rho*e/Ln*expo_age*1E+06));
    N10ms_pb=P10ms_z/(1/tau_10/1E+06+rho*e/Lms)*(1-exp(-expo_age/tau_10-rho*e/Lms*expo_age*1E+06));
    N10mf_pb=P10mf_z/(1/tau_10/1E+06+rho*e/Lmf)*(1-exp(-expo_age/tau_10-rho*e/Lmf*expo_age*1E+06));
    N10pb=N10n_pb+N10ms_pb+N10mf_pb;
    
    % Monte-Carlo for the error of post-burial 10Be
    cache=zeros(1,simulation_times);
    for i=1:simulation_times
        rand_P10n_z=normrnd(P10n_z,sigma_P10n_z);
        rand_P10ms_z=normrnd(P10ms_z,sigma_P10ms_z);
        rand_P10mf_z=normrnd(P10mf_z,sigma_P10mf_z);
        rand_tau_10=normrnd(tau_10,sigma_tau_10);
        rand_expo_age=expo_age_est(randi([1,n]));
        if option.flag2==1  % use erosion rate
            rand_e=e_est(randi([1,n2]));
            rand_N10n_pb=rand_P10n_z/(1/rand_tau_10/1E+06+rho*rand_e/Ln)*(1-exp(-rand_expo_age/rand_tau_10-rho*rand_e/Ln*rand_expo_age*1E+06));
            rand_N10ms_pb=rand_P10ms_z/(1/rand_tau_10/1E+06+rho*rand_e/Lms)*(1-exp(-rand_expo_age/rand_tau_10-rho*rand_e/Lms*rand_expo_age*1E+06));
            rand_N10mf_pb=rand_P10mf_z/(1/rand_tau_10/1E+06+rho*rand_e/Lmf)*(1-exp(-rand_expo_age/rand_tau_10-rho*rand_e/Lmf*rand_expo_age*1E+06));
        elseif option.flag2==0  % no erosion rate
            rand_N10n_pb=rand_P10n_z/(1/rand_tau_10/1E+06+rho*e/Ln)*(1-exp(-rand_expo_age/rand_tau_10-rho*e/Ln*rand_expo_age*1E+06));
            rand_N10ms_pb=rand_P10ms_z/(1/rand_tau_10/1E+06+rho*e/Lms)*(1-exp(-rand_expo_age/rand_tau_10-rho*e/Lms*rand_expo_age*1E+06));
            rand_N10mf_pb=rand_P10mf_z/(1/rand_tau_10/1E+06+rho*e/Lmf)*(1-exp(-rand_expo_age/rand_tau_10-rho*e/Lmf*rand_expo_age*1E+06));
        end
        rand_N10pb=rand_N10n_pb+rand_N10ms_pb+rand_N10mf_pb;
        cache(i)=(rand_N10pb);
    end
    sigma_N10pb=std(cache);

    % post-burial 26Al (eqn6 in Lal, 1991)
    N26n_pb=P26n_z/(1/tau_26/1E+06+rho*e/Ln)*(1-exp(-expo_age/tau_26-rho*e/Ln*expo_age*1E+06));
    N26ms_pb=P26ms_z/(1/tau_26/1E+06+rho*e/Lms)*(1-exp(-expo_age/tau_26-rho*e/Lms*expo_age*1E+06));
    N26mf_pb=P26mf_z/(1/tau_26/1E+06+rho*e/Lmf)*(1-exp(-expo_age/tau_26-rho*e/Lmf*expo_age*1E+06));
    N26pb=N26n_pb+N26ms_pb+N26mf_pb;

    % Monte-Carlo for the error of post-burial 26Al
    cache=zeros(1,simulation_times);
    for i=1:simulation_times
        rand_P26n_z=normrnd(P26n_z,sigma_P26n_z);
        rand_P26ms_z=normrnd(P26ms_z,sigma_P26ms_z);
        rand_P26mf_z=normrnd(P26mf_z,sigma_P26mf_z);
        rand_tau_26=normrnd(tau_26,sigma_tau_26);
        rand_expo_age=expo_age_est(randi([1,n]));
        if option.flag2==1  % use erosion rate
            rand_e=e_est(randi([1,n2]));
            rand_N26n_pb=rand_P26n_z/(1/rand_tau_26/1E+06+rho*rand_e/Ln)*(1-exp(-rand_expo_age/rand_tau_26-rho*rand_e/Ln*rand_expo_age*1E+06));
            rand_N26ms_pb=rand_P26ms_z/(1/rand_tau_26/1E+06+rho*rand_e/Lms)*(1-exp(-rand_expo_age/rand_tau_26-rho*rand_e/Lms*rand_expo_age*1E+06));
            rand_N26mf_pb=rand_P26mf_z/(1/rand_tau_26/1E+06+rho*rand_e/Lmf)*(1-exp(-rand_expo_age/rand_tau_26-rho*rand_e/Lmf*rand_expo_age*1E+06));
        elseif option.flag2==0  % no erosion rate
            rand_N26n_pb=rand_P26n_z/(1/rand_tau_26/1E+06+rho*e/Ln)*(1-exp(-rand_expo_age/rand_tau_26-rho*e/Ln*rand_expo_age*1E+06));
            rand_N26ms_pb=rand_P26ms_z/(1/rand_tau_26/1E+06+rho*e/Lms)*(1-exp(-rand_expo_age/rand_tau_26-rho*e/Lms*rand_expo_age*1E+06));
            rand_N26mf_pb=rand_P26mf_z/(1/rand_tau_26/1E+06+rho*e/Lmf)*(1-exp(-rand_expo_age/rand_tau_26-rho*e/Lmf*rand_expo_age*1E+06));
        end
        rand_N26pb=rand_N26n_pb+rand_N26ms_pb+rand_N26mf_pb;
        cache(i)=(rand_N26pb);
    end
    sigma_N26pb=std(cache);
end
