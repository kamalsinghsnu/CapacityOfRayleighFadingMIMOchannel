% The Fig. 9 in [1] is verified/checked for errors/discrepancies
%
% Written by Kamal Singh, kamal.singh@snu.edu.in
%==========================================================================


clear all;
clc;            



global m n;


% Choose parameters "m", "n" and SNR (in dB)
%------------------------------------------------------
m = 2; % same as fixed in Fig. 9 in [1] 


n_range = [2;3;4;5;8]; % same as fixed in Fig. 9 in [1] 


len_n = length(n_range);


SNR_dB =  -10:5:10;  % in dB. This is used as SNR in the Fig. 9 plot in [1] 


len = length(SNR_dB);


Outage_prob_bound = zeros(len,1);  % Numerical Capacity estimates (with CSIT) 


for j = 1:len_n

    n = n_range(j);
    
    fid1 = fopen(sprintf('Outage probability for (Fig. 9 in Jayaweera paper) for %d x %d Rayleigh channel (Full CSIT).txt',m,n),'wt');

for i = 1:1:len
            
                
                
                SNR = 10^(SNR_dB(i) / 10); % SNR in linear scale
            

                


                
                % Computing the Waterfilling threshold lamda0 
                overall = @(x) integral(@(y) (1./x - 1./y) .*  pdf_lamda(y),x,inf) - SNR;
                lamda0 = fzero(overall,[1e-100,1000]);

             
                % Outage Probability Upperbound computation: Eq. (64) in [1]
                Outage_prob_bound(i) =  (gamma(n+m-1) -  gamma(n+m-1) .* gammainc(lamda0,n+m-1,'upper')) ./(gamma(n)*gamma(m));
                
                fprintf(fid1,'(%.15f,%.15f)\n',SNR_dB(i),Outage_prob_bound(i));
                
                               
                
end

fclose(fid1);   

plot(SNR_dB,Outage_prob_bound,'r');

hold on;

end