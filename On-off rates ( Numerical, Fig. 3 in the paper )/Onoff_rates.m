% This MATLAB script generates a numerical estimate of the Ergodic
% rates achievable with On-off power scheme on the Best-eigenmode 
% of the coherent Rayleigh fading MIMO channel "with CSIT" of 
% dimension Nr x Nt.
% 
%
% At the end of this routine, the On-off rate values are stored in 
% the variable Onoff_rates -- see below
%
% For comparison purpose, the script also creates a numerical estimate 
% of the Ergodic Capacity of the coherent Rayleigh fading MIMO channel 
% "with CSIT" of dimension Nr x Nt: the values are stored in the
% variable Capacity_CSIT -- see below
%
% Written by Kamal Singh, kamal.singh@snu.edu.in
%==========================================================================

clear all;
clc;

global m n lamda0;

Nr = 4; % No. of receive antennas

Nt = 10; % No. of transmit antennas

m = min(Nr,Nt);

n = max(Nr,Nt);

SNR_dB =  0:-5:-40;  % in dB. This is used as SNR in the plots in the paper 

len = length(SNR_dB);





Onoff_rates = zeros(len,1); % storing On-off rates (numerical estimates)

Capacity_CSIT = zeros(len,1); % Waterfilling rates = Capacity (numerical estimates)








for i = 1:1:len
            
            
            

                SNR = 10^(SNR_dB(i) / 10); % linear scale

                
                % Computing the Waterfilling threshold lamda0
                overall = @(x) integral(@(y) (1./x - 1./y) .*  pdf_lamda(y),x,inf) - SNR;
                lamda0 = fzero(overall,[1e-100,1000]);
            
        
                % On-off power control parameter "tau"
                tau = lamda0;  % On-off cutoff threshold taken same as Waterfilling cutoff threshold
                
                P0 = m .* SNR ./ (1 - CDF_max(tau)); % Tall & Alouini paper (SEction IV)
      
                
                % Ergodic-Rates for the with on-off power control on the
                % "best-user" are next computed numerically
                Onoff_rate_integral = @(y) (1-CDF_max(y)) .* P0 ./ (1 + P0 .* y);
                Onoff_rates(i) = log2(2.7182818) .* ((1-CDF_max(tau)) .* log(1 + P0 .* tau) + integral(Onoff_rate_integral,tau,1000,'ArrayValued',true));
                 
                
                % Ergodic capacity with Waterfilling scheme (for comparison)
                Capacity_integral = @(v) log2(v/lamda0) .* pdf_lamda(v);
                Capacity_CSIT(i) = m * integral(Capacity_integral,lamda0,inf);
                
                
                [SNR_dB'     Capacity_CSIT        Onoff_rates]


end
   
      

     
plot(SNR_dB,Onoff_rates,'g');
hold on;

plot(SNR_dB,Capacity_CSIT,'r');
hold on;
    
    
