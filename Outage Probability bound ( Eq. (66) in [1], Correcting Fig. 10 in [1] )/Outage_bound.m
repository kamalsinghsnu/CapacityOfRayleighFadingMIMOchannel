% The Fig. 10 in [1] is verified/checked for errors/discrepancies
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


Outage_prob_bound = zeros(len,1);  % Numerical Capacity estimates (with CSIT) 


for j = 1:len_n

    n = n_range(j);
    
    m = n; %  In Fig. 10 in [1], this setting is selected
    
    fid1 = fopen(sprintf('Outage probability for (Fig. 10 in Jayaweera paper) for %d x %d Rayleigh channel (Full CSIT).txt',m,n),'wt');

for i = 1:1:len
            
                
                
                SNR = 10^(SNR_dB(i) / 10); % SNR in linear scale
            

                


                
                % Computing the Waterfilling threshold lamda0 
                overall = @(x) integral(@(y) (1./x - 1./y) .*  pdf_lamda(y),x,inf) - SNR;
                lamda0 = fzero(overall,[1e-100,1000]);

             
                % Outage Probability Upperbound computation: Eq. (64) in [1]
                % gammainc function in MATLAB is a scaled (divided by gamma(a)) gamma function,
                % therefore, to perform Eq. (64), we remove this scaling by
                % multiplying by gamma(a). 
                p1 = (gamma(n+m-1) -  gamma(n+m-1) .* gammainc(lamda0,n+m-1,'upper')) ./(gamma(n)*gamma(m)); %Eq. (64) upper bound 1
                p2 = 1 -exp(-m*lamda0); %Eq. (66) upper bound 2
                
                Outage_prob_bound(i) = min(p1,p2);
                
                fprintf(fid1,'(%d,%.30f)\n',SNR_dB(i),Outage_prob_bound(i));
                
                
                %NOTE: For a small subset of SNR and m, n parameters, the
                %code is generating zero values becau
                
                               
                
end

fclose(fid1); 

plot(SNR_dB,Outage_prob_bound,'r');

hold on;

end