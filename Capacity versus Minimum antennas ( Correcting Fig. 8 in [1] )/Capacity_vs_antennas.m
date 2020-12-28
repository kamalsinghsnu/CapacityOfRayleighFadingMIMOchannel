% The Fig. 8 in [1] is verified/checked for errors/discrepancies
%
% Written by Kamal Singh, kamal.singh@snu.edu.in
%==========================================================================

clear all;
clc;            

global m n;

% Choose parameters "m", "n" and SNR (in dB)
%------------------------------------------------------
n = 18; % same as fixed in Fig. 8 in [1] 


m_range = [2:1:8]; % same as fixed in Fig. 8 in [1] 


len_m = length(m_range);


SNR_dB =  0:5:20;  % in dB. This is used as SNR in the Fig. 8 plot in [1] 


len = length(SNR_dB);


Capacity_CSIT = zeros(len_m,1);  % Numerical Capacity estimates (with CSIT) 



for i = 1:len

    SNR = 10^(SNR_dB(i) / 10); % SNR in linear scale
    
    fid1 = fopen(sprintf('Capacity versus minimum antennas for SNR = %d dB (Fig. 8 in Jayaweera paper) for n = %d Rayleigh channel (Full CSIT).txt',SNR_dB(i),n),'wt');

for j = 1:1:len_m
            
                m = m_range(j);
            

                


                
                % Computing the Waterfilling threshold lamda0 
                overall = @(x) integral(@(y) (1./x - 1./y) .*  pdf_lamda(y),x,inf) - SNR;
                lamda0 = fzero(overall,[1e-100,1000]);

             
                % Ergodic capacity With CSIT
                Capacity_integral = @(v) log2(v/lamda0) .* pdf_lamda(v);
                Capacity_CSIT(j) = m * integral(Capacity_integral,lamda0,inf);
                
                fprintf(fid1,'(%d,%.6f)\n',m,Capacity_CSIT(j));
                
                               
                
end

fclose(fid1);

plot(m_range,Capacity_CSIT,'r');

hold on;

end