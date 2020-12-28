% The Fig. 3 in the paper [1] is being checked for any errors if present
%=====================================================================================
% This MATLAB script generates a numerical estimate of the Ergodic
% rates achievable with Maximal Gain transmission with Waterfilling 
% power scheme the 1 x Nt (i.e. SINGLE ANTENNA RECEIVER SYSTEMS) 
% Rayleigh fading MIMO channel "with CSIT".
% 
%
% At the end of this routine, the "maximal gain rates" values are stored 
% in the variable Maximal_Gain_rates -- see below
%
% Written by Kamal Singh, kamal.singh@snu.edu.in
%==========================================================================

clear all;
clc;

global m n lamda0;

Nr = 1; % No. of receive antennas

Nt = 8; % No. of transmit antennas

m = min(Nr,Nt);

n = max(Nr,Nt);

SNR_dB =  -10:5:15;  % in dB. This is used as SNR in the plots in the paper 

len = length(SNR_dB);





Maximal_Gain_rates = zeros(len,1); % storing maximal gain transmission=rates (numerical estimates)


% PDF of the largest of all the magnitude-squared fading coefficients
pdf_lamda_max = @(x) n .* (1 - exp(-x)).^(n-1) .* exp(-x);



for i = 1:1:len
            
            
            

                SNR = 10^(SNR_dB(i) / 10); % linear scale

                
                % Computing the Waterfilling threshold lamda0
                overall = @(x) integral(@(y) (1./x - 1./y) .*  pdf_lamda_max(y),x,inf) - SNR;
                lamda0 = fzero(overall,[1e-100,1000]);
            
        
                % Maximal Gain Transmission-rates with Waterfilling scheme (for comparison)
                Capacity_integral = @(v) log2(v/lamda0) .* pdf_lamda_max(v);
                Maximal_Gain_rates(i) = integral(Capacity_integral,lamda0,inf);
                
                
                [SNR_dB'     Maximal_Gain_rates]


end
   
      

     
plot(SNR_dB,Maximal_Gain_rates,'b');
hold on;

    
    
