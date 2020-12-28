
% The Fig. 1 and Fig. 4 in the paper [1] is being checked for any errors if present
%=====================================================================================
% This MATLAB script generates the optimal cut-off of the Waterfilling power 
% scheme versus SNR, which is optimal power control for the coherent Rayleigh 
% fading MIMO channel "with CSIT".
% 
%
% At the end of this routine, the resulting cut-off values are stored in 
% the variable Optimal_cutoff -- see below
%
% 
% Written by Kamal Singh, kamal.singh@snu.edu.in
%=====================================================================================



clear all;
clc;            

global m n;

% Choose parameters "m" and "n"

m = 4; % For Fig.1, set m = 1 (rank 1 (one) channel)

n_vector = 4:2:10;

len_n = length(n_vector);

% Choose SNR (in dB) for the desired Figure (in the paper [1])

SNR_dB =  -10:5:40;  % in dB. This is used as SNR in the plots in the paper 

len = length(SNR_dB);



Optimal_cutoff = zeros(len,1);  % Waterfilling cutoff estimates 



for h = 1:1:len_n 

    n = n_vector(h);
    
    %fid = fopen(sprintf('Gamma_0 readings for %d x %d Rayleigh channel (Full CSIT).txt',m,n),'wt');

for i = 1:1:len
            
            
            

                SNR = 10^(SNR_dB(i) / 10); % SNR in linear scale


                
                % Computing the threshold/cut-off in Jayaweera's
                % paper [1] (Fig. 4 and Fig. 1) using the pdf of the unordered Eigenvalue. 
                overall = @(x) integral(@(y) (1./x - 1./(y .* SNR )) .*  pdf_lamda(y),x./SNR,inf) - 1;
                Optimal_cutoff(i) = fzero(overall,[1e-100,1000]);

                
                % Saving the readings of Cut-off versus SNR (in dB)
                %fprintf(fid,'(%.6f,%.6f)\n',SNR_dB(i),Optimal_cutoff(i));
                
                
                [SNR_dB'   Optimal_cutoff]




end

%fclose(fid);

plot(SNR_dB,Optimal_cutoff);

hold on;

end