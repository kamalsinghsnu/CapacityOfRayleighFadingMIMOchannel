% Plots in the Fig. 9 in the submitted Manuscript are generated
%==========================================================================
% We compute the actual Outage Probability (P_out) and compare with the
% corrected Upper bound (Eq. (66) in [1])
% 
% Written by Kamal Singh, kamal.singh@snu.edu.in
%==========================================================================


clear all;
clc;            


global m n;


% Choose parameters "m", "n" and SNR (in dB)
%------------------------------------------------------


n_range = [2;3;4;6;8]; % same range as fixed in Fig. 10 in [1] 

len_n = length(n_range);

SNR_dB =  [-2;0;5;10];  % in dB. This is used as SNR in the Fig. 10 plot in [1] 




len = length(SNR_dB);




Outage_prob_bound = zeros(len,1);  %  Eq. (66) in [1]
Outage_prob_actual = zeros(len,1);  %  Based on CDF of the largest Eigen-value


for j = 1:len_n

    n = n_range(j);
    
    m = n; %  In Fig. 10 in [1], this setting is selected

for i = 1:1:len
            
                
                
                SNR = 10^(SNR_dB(i) / 10); % SNR in linear scale
            

                


                
                % Computing the Waterfilling threshold lamda0 
                overall = @(x) integral(@(y) (1./x - 1./y) .*  pdf_lamda(y),x,inf) - SNR;
                lamda0 = fzero(overall,[1e-100,1000]);

             
                % Outage Probability Upperbound computation: Eq. (64) in [1]
                p1 =  (gamma(n+m-1) -  gamma(n+m-1) .* gammainc(lamda0,n+m-1,'upper')) ./(gamma(n)*gamma(m)); %Eq. (64) upper bound 1
                                
                % Eq. (65) in [1] needed for Eq. (66)
                p2 = 1 -exp(-m*lamda0); 
                
                
                % Outage Probability Upperbound computation: Eq. (66) in [1]
                Outage_prob_bound(i) = min(p1,p2); 
                
                
                Outage_prob_actual(i) = CDF_max(lamda0);
                
                
                
                
                
                
                
                               
                
end

% RED for the bound (Eq. (66))
plot(SNR_dB,Outage_prob_bound,'r');

hold on;

% BLUE for the actual
plot(SNR_dB,Outage_prob_actual,'b');

hold on;

end