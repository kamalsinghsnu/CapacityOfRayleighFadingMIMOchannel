% This MATLAB script generates a Monte Carlo estimate of the Ergodic
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

global m n;


samples = 10^6; % No. of samples for Monte-carlo simulation

Nr = 4; % No. of receive antennas

Nt = 10; % No. of transmit antennas

m = min(Nr,Nt);

n = max(Nr,Nt);

SNR_dB =  0:-5:-40;  % in dB. This is used as SNR in the plots in the paper 

len = length(SNR_dB);






Onoff_rates = zeros(len,1); % storing On-off rates (Monte Carlo estimates)

Capacity_CSIT = zeros(len,1); % Waterfilling rates = Capacity (numerical estimates)









for i = 1:1:len
            
            
            

                SNR = 10^(SNR_dB(i) / 10); % SNR in linear scale

                
                % Computing the Waterfilling threshold lamda0
                overall = @(x) integral(@(y) (1./x - 1./y) .*  pdf_lamda(y),x,inf) - SNR;
                lamda0 = fzero(overall,[1e-100,1000]);
                
                
                % Waterfilling rates (=capacity) computed numerically for
                % for comparison purposes
                Capacity_integral = @(v) log2(v/lamda0) .* pdf_lamda(v);
                Capacity_CSIT(i) = m * integral(Capacity_integral,lamda0,inf);
                
                
                
                
                % We need to use the Largest Eigenvalue CDF as will be
                % shown next
                   
                
                
                % CDF of the Largest EigenValue of the Wishart Matrix is
                % obtained below
                
                funcky = @(a,b,x)  gammainc(x,n-m+a+b-1,'lower') .* gamma(n-m+a+b-1);
                
                divider = 1;
                
                for kk = 1:1:m
                    
                    divider = divider * gamma(n-kk+1)* gamma(m-kk+1);
                    
                end
                
                numerator_matrix = zeros(m,m);
                
                cutoff = lamda0;
                
                for ii = 1:1:m
                    
                    for jj = 1:1:m
                        
                        
                        numerator_matrix(ii,jj) = funcky(ii,jj,cutoff);
                        
                        
                    end
                    
                end
                
                
                CDF_tau = det(numerator_matrix) ./ divider;
                
                
                
                
                
                
                
                
                
                % The On-off Power level
                P0 = SNR * m ./ (1 - CDF_tau); 
                
   


               
                % Monte-carlo simulation starts from here onwards 
                Temp = 0;
                
                
                
              
                for k=1:1:samples
                    
                    
                                H = (randn(Nr,Nt)+1i*randn(Nr,Nt))./sqrt(2);
                                
                                [U,S,V] = svd(H);
                                
                                
                                lamda_max = max(diag(S));
                                
                                                                                              
                                if(lamda_max^2>= cutoff)
                                    Temp = Temp  + log2(1 + P0 * lamda_max^2);
                                end
                                
                end
                
                
                Onoff_rates(i) = Temp ./ samples;
                
                %Simulation ends here...
                
                
                
                [SNR_dB'     Capacity_CSIT     Onoff_rates]
                
                


end
    
    
    

    
    