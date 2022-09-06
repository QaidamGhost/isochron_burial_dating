function [Pn_z,sigma_Pn_z,Pms_z,sigma_Pms_z,Pmf_z,sigma_Pmf_z] = production_rate(measured_lat,measured_elv,shielding_factor,z,rho,nuclide)

%% Calculate 10Be or 26Al Production Rates on the surface or at a given depth from Spallation and Muons.
% This script is modified after a MATLAB script "CNP.m" which was developed
% by Chia-Yu Chen, Richard Ott, Erica Erlanger, Maarten Lupker, and Yanyan
% Wang (Lupker et al., 2012). For a surface production rate calculation,
% just set the "z" to zero and "rho" to an arbitrary value.
% Note that the Spallation-induced production rates ratio here is 6.61.

%% Arguments:
% measured_lat: mearsured latitude of the sample (unitless; scalar)
% measured_elv: mearsured elevation of the sample (m; scalar)
% shielding_factor: (unitless; scalar)
% z: measured depth of the sample (cm; scalar)
% rho: density of the overburden (g/cm^3; scalar)
% nuclide: value 10 for 10Be calculation and 26 for 26Al (unitless; scalar)

%% Output:
% Pn_z: production rate from spallation at the depth (atom/g/yr; scalar)
% sigma_Pn_z: 1 sigma error of Pn_z (atom/g/yr; scalar)
% Pms_z: production rate from slow muons at the depth (atom/g/yr; scalar)
% sigma_Pms_z: 1 sigma error of Pms_z (atom/g/yr; scalar)
% Pmf_z: production rate from fast muons at the depth (atom/g/yr; scalar)
% sigma_Pmf_z: 1 sigma error of Pmf_z (atom/g/yr; scalar)

    if nuclide~=10&&nuclide~=26 % wrong input of the "nuclide"
        Pn_z=0;
        sigma_Pn_z=0;
        Pms_z=0;
        sigma_Pms_z=0;
        Pmf_z=0;
        sigma_Pmf_z=0;
        return;
    end

    %% Scaling Equation Constants by Latitude (Stone,2000 after Lal, 1991)
    load consts.mat Lat a b c d e m;
    % Interpolation
    A = interp1(Lat,a,measured_lat);
    B = interp1(Lat,b,measured_lat);
    C = interp1(Lat,c,measured_lat);
    D = interp1(Lat,d,measured_lat);
    E = interp1(Lat,e,measured_lat);
    M = interp1(Lat,m,measured_lat);

    %% Perform air pressure correction for elevation (Eq. 1 from Stone, 2000)to account for diminished cosmic ray attenuation with decreased atmospheric pressue 
    pres=1013.25*exp(((-0.03417)/6.5e-3)*(log(288.15)-log(288.15-(6.5e-3*measured_elv)))); 

    %% Some values for 10Be and 26Al production rates at SLHL in Braucher et al., 2011
    load consts.mat Pn10_SLHL sigma_Pn10_SLHL Pms10_SLHL sigma_Pms10_SLHL Pmf10_SLHL sigma_Pmf10_SLHL Pn26_SLHL sigma_Pn26_SLHL Pms26_SLHL sigma_Pms26_SLHL Pmf26_SLHL sigma_Pmf26_SLHL;
    if nuclide==10
        % Spallation-induced Production Rate of 10Be (at/g/yr) assuming 'St" scaling framework (Braucher et al., 2011)
        Pn_SLHL = Pn10_SLHL; 
        sigma_Pn_SLHL = sigma_Pn10_SLHL;
        % Slow Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al., 2011)
        Pms_SLHL = Pms10_SLHL; 
        sigma_Pms_SLHL = sigma_Pms10_SLHL;
        % Fast Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al., 2011)
        Pmf_SLHL = Pmf10_SLHL; 
        sigma_Pmf_SLHL=sigma_Pmf10_SLHL;
    elseif nuclide==26
        % Spallation-induced Production Rate of 26Al (at/g/yr) assuming 'St" scaling framework (Braucher et al., 2011)
        Pn_SLHL = Pn26_SLHL; 
        sigma_Pn_SLHL = sigma_Pn26_SLHL; 
        % Slow Muon-induced Production Rate of 26Al (at/g/yr) scaled to sea level (Braucher et al., 2011)
        Pms_SLHL = Pms26_SLHL; 
        sigma_Pms_SLHL = sigma_Pms26_SLHL; 
        % Fast Muon-induced Production Rate of 26Al (at/g/yr) scaled to sea level (Braucher et al., 2011)
        Pmf_SLHL = Pmf26_SLHL; 
        sigma_Pmf_SLHL=sigma_Pmf26_SLHL;
    end

    %% Calculate 10Be or 26Al Production Rates on the surface from Spallation and Muons
    % 10Be or 26Al production from Spallation, assuming neutron attenuation length in air of 150 g/cm2
    Pn= Pn_SLHL.*(A+B.*exp(-pres/150)+C.*pres+D.*pres.^2+E.*pres.^3)*shielding_factor; % Stone, 2000
    sigma_Pn= sigma_Pn_SLHL.*(A+B.*exp(-pres/150)+C.*pres+D.*pres.^2+E.*pres.^3)*shielding_factor;
    % 10Be or 26Al production from Slow Muons, assuming (1) sea level pressure of 1013.25 mbar and 
    %(2) muon attentuation length in air of 260 g/cm2 (Braucher et al., 2011)
    Pms=Pms_SLHL.*exp((1013.25-pres)/260)*shielding_factor; 
    sigma_Pms= sigma_Pms_SLHL.*exp((1013.25-pres)/260)*shielding_factor; 
    % 10Be or 26Al production from Fast Muons, assuming (1) sea level pressure of 1013.25 mbar and 
    %(2) muon attentuation length in air of 510 g/cm2 (Braucher et al., 2011)
    Pmf=Pmf_SLHL.*exp((1013.25-pres)/510)*shielding_factor; 
    sigma_Pmf=sigma_Pmf_SLHL.*exp((1013.25-pres)/510)*shielding_factor; 

    %% Attenuation lengths of neutrons, slow muons, and fast muons, respectively. 
    % Unit: g/cm^2. Braucher et al., 2011.
    load consts.mat Ln Lms Lmf;
    
    %% Calculate 10Be or 26Al Production Rates at a given depth from Spallation and Muons
    Pn_z= Pn*exp(-z*rho/Ln);
    sigma_Pn_z= sigma_Pn*exp(-z*rho/Ln);
    Pms_z= Pms*exp(-z*rho/Lms);
    sigma_Pms_z= sigma_Pms*exp(-z*rho/Lms);
    Pmf_z= Pmf*exp(-z*rho/Lmf);
    sigma_Pmf_z= sigma_Pmf*exp(-z*rho/Lmf);
end
