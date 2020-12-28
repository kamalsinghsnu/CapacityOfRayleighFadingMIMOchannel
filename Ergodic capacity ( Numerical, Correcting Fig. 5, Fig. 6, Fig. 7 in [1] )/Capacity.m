% This MATLAB script generates a numerical estimate of the Ergodic
% Capacity of the coherent Rayleigh fading MIMO channel "with CSIT" of 
% dimension Nr x Nt.
% 
%
% At the end of this routine, the resulting capacity values are stored in 
% the variable Capacity_CSIT -- see below
%
% For comparison purpose, the script also creates a numerical estimate 
% of the Ergodic Capacity of the coherent Rayleigh fading MIMO channel 
% "without CSIT" of dimension Nr x Nt: the values are stored in the
% variable Capacity_noCSIT -- see below
%
% Written by Kamal Singh, kamal.singh@snu.edu.in
%==========================================================================

clear all;
clc;            

global m n;


Nr = 4; % No. of receive antennas

Nt = 10; % No. of transmit antennas

m = min(Nr,Nt);

n = max(Nr,Nt);

SNR_dB =  -10:5:15;  % in dB. This is used as SNR in the plots in the paper 

len = length(SNR_dB);



Capacity_CSIT = zeros(len,1);  % Numerical Capacity estimates (with CSIT) 
Capacity_noCSIT = zeros(len,1); % Numerical Capacity estimates (without CSIT) 




for i = 1:1:len
            
            
            

                SNR = 10^(SNR_dB(i) / 10); % SNR in linear scale


                
                % Computing the Waterfilling threshold lamda0 
                overall = @(x) integral(@(y) (1./x - 1./y) .*  pdf_lamda(y),x,inf) - SNR;
                lamda0 = fzero(overall,[1e-100,1000]);

             
                % Ergodic capacity With CSIT
                Capacity_integral = @(v) log2(v/lamda0) .* pdf_lamda(v);
                Capacity_CSIT(i) = m * integral(Capacity_integral,lamda0,inf);
                
                % Ergodic capacity Without CSIT
                Capacity_noCSIT(i) =  m * integral(@(v) log2(1 + v .* m .* SNR./ Nt) .* pdf_lamda(v),0,inf);
                
                
                [SNR_dB'   Capacity_CSIT    Capacity_noCSIT]




end

      
plot(SNR_dB,Capacity_CSIT,'r');

hold on;

plot(SNR_dB,Capacity_noCSIT,'b');

hold on;