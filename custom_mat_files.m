function custom_mat_files(erosion_rate,expo_age,sigma_expo_age)

%% Generate "e.mat" and "expo_age.mat" if you want to calculate a maximum isochron burial age and already know the value and absolate 1 sigma error of erosion rate and exposure age.
% Note that the input exposure age should be Gaussian distribution, or you
% should use modified Hidy's scripts to generate the "expo_age.mat" file.

expo_age_est=zeros(1E4,1);
for i=1:1E4
  expo_age_est(i)=normrnd(expo_age,sigma_expo_age);
end
save expo_age.mat expo_age expo_age_est;
e=erosion_rate;
save e.mat e;
end