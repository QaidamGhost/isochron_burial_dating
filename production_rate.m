function [Pn_z,sigma_Pn_z,Pms_z,sigma_Pms_z,Pmf_z,sigma_Pmf_z] = production_rate(measured_lat,measured_elv,shielding_factor,z,rho,nuclide)

%% Calculate 10Be or 26Al Production Rates on the surface or at a given depth from Spallation and Muons.
% This script is modified after a MATLAB script "CNP.m" which was developed
% by Chia-Yu Chen, Richard Ott, Erica Erlanger, Maarten Lupker, and Yanyan
% Wang (Lupker et al., 2012). For a surface production rate calculation,
% just set the "z" to zero and "rho" to an arbitrary value.

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

    %% Scaling Equation Constants by Latitude (Stone,2000)
    Lat = [0,10,20,30,40,50,60,90];
    a = [31.8518,34.3699,40.3153,42.0983,56.7733,69.0720,71.8733,71.8733];
    b = [250.3193,258.4759,308.9894,512.6857,649.1343,832.4566,863.1927,863.1927];
    c = [-0.083393,-0.089807,-0.106248,-0.120551,-0.160859,-0.199252,-0.207069,-0.207069];
    d = [7.4260e-5,7.9457e-5,9.4508e-5, 1.1752e-4,1.5463e-4,1.9391e-4,2.0127e-4,2.0127e-4];
    e = [-2.2397e-8,-2.3697e-8,-2.8234e-8,-3.8809e-8,-5.0330e-8,-6.3653e-8,-6.6043e-8,-6.6043e-8];
    m = [0.587,0.600,0.678,0.833,0.933,1.000,1.000,1.000];
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
    if nuclide==10
        % Spallation-induced Production Rate of 10Be (at/g/yr) assuming 'St" scaling framework (Braucher et al., 2011)
        Pn_SLHL = 4.49; 
        sigma_Pn_SLHL = 0.30;
        % Slow Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al., 2011)
        Pms_SLHL = 0.012; 
        sigma_Pms_SLHL = 0.012;
        % Fast Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al., 2011)
        Pmf_SLHL = 0.039; 
        sigma_Pmf_SLHL=0.004;
    elseif nuclide==26
        % Spallation-induced Production Rate of 26Al (at/g/yr) assuming 'St" scaling framework (Braucher et al., 2011)
        Pn_SLHL = 29.68; 
        sigma_Pn_SLHL = 1.98; 
        % Slow Muon-induced Production Rate of 26Al (at/g/yr) scaled to sea level (Braucher et al., 2011)
        Pms_SLHL = 0.84; 
        sigma_Pms_SLHL = 0.17; 
        % Fast Muon-induced Production Rate of 26Al (at/g/yr) scaled to sea level (Braucher et al., 2011)
        Pmf_SLHL = 0.081; 
        sigma_Pmf_SLHL=0.051;
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
    Ln=160;
    Lms=1500;
    Lmf=4320;
    
    %% Calculate 10Be or 26Al Production Rates at a given depth from Spallation and Muons
    Pn_z= Pn*exp(-z*rho/Ln);
    sigma_Pn_z= sigma_Pn*exp(-z*rho/Ln);
    Pms_z= Pms*exp(-z*rho/Lms);
    sigma_Pms_z= sigma_Pms*exp(-z*rho/Lms);
    Pmf_z= Pmf*exp(-z*rho/Lmf);
    sigma_Pmf_z= sigma_Pmf*exp(-z*rho/Lmf);
end
