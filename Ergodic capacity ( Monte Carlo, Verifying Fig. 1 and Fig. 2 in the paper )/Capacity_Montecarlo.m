% This MATLAB script generates a Monte Carlo estimate of the Ergodic
% Capacity of the coherent Rayleigh fading MIMO channel "with CSIT" of 
% dimension Nr x Nt.
% 
%
% At the end of this routine, the resulting capacity values are stored in 
% the variable Capacity_CSIT -- see below
%
% For comparison purpose, the script also creates a Monte Carlo estimate 
% of the Ergodic Capacity of the coherent Rayleigh fading MIMO channel 
% "without CSIT" of dimension Nr x Nt: the values are stored in the
% variable Capacity_noCSIT -- see below
%
% Written by Kamal Singh, kamal.singh@snu.edu.in
%==========================================================================

clear all;
clc;

global m n;

samples = 10^6; % No. of samples for Monte-carlo simulation

Nr = 1; % No. of receive antennas

Nt = 10; % No. of transmit antennas

m = min(Nr,Nt);

n = max(Nr,Nt);

SNR_dB =  -15:5:10;  % in dB. This is used as SNR in the plots in the paper 


len = length(SNR_dB);



Capacity_CSIT = zeros(len,1);  % Monte-Carlo Capacity estimates (with CSIT) 
Capacity_noCSIT = zeros(len,1); % Monte-Carlo Capacity estimates (without CSIT) 



for i = 1:1:len
            
            
            

                SNR = 10^(SNR_dB(i) / 10); % SNR in linear scale



                % Computing the Waterfilling threshold lamda0 
                overall = @(x) integral(@(y) (1./x - 1./y) .*  pdf_lamda(y),x,inf) - SNR;
                lamda0 = fzero(overall,[1e-100,1000]);



              
                % Monte-carlo simulation starts from here onwards 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                Temp = 0; % Final value will store capacity value with CSIT
                                
                Temp1 = 0; % Final value will store capacity value without CSIT
                
                for k=1:1:samples
                    
                                %%%%%%%%%%%%%%%%%%%%%%Channel 1 realizations%%%%%%%%%%
                                H_1 = (randn(Nr,Nt)+1i*randn(Nr,Nt))./sqrt(2);
                                
                                [U,S,V] = svd(H_1);
                                perm = randperm(m);
                                Delta = zeros(Nr,Nt);
                                %%%%%%%%%%%%%%% Delta1 permuted %%%%%%%%%%%%%%%%%%%%
                                for z = 1:1:m,
                                   Delta(z,z) = diag(S(perm(z),perm(z)));
                                end
                                %%%%%%%%%%%%%%  U1 permuted next %%%%%%%%%%%%%%%%%%%
                                U1 = zeros(Nr,Nr);
                                for z = 1:1:m,
                                    U1(:,z) = U(:,perm(z));
                                end
                                if(m~=Nr)
                                    for z = m+1:1:Nr,
                                        U1(:,z) = U(:,z);
                                    end
                                end


                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%% With CSIT: Waterfilling power control
                                Q = zeros(Nt,Nt);
                                for kk = 1:1:m
                                    if(Delta(kk,kk) >=  sqrt(lamda0))
                                        Q(kk,kk) = max(((1./lamda0) - (1./(Delta(kk,kk).^2))),0);
                                    else
                                        Q(kk,kk) = 0;
                                    end
                                end
                                
                                X_Q = eye(Nr) + Delta * Q * Delta';
                                det_X_Q = real(det(X_Q));
                                Temp = Temp  + log2(det_X_Q);
                                
                                
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
                                %%%%%%%%%%%%%% Without CSIT power control
                                Q_noCSIT = (m * SNR/Nt) .* diag(ones(Nt,1));
                                
                
                                
                                X_Q_noCSIT = eye(Nr) + U1 * Delta * Q_noCSIT * Delta' * U1' ;
                                det_X_Q_noCSIT = real(det(X_Q_noCSIT));
                                Temp1 = Temp1  + log2(det_X_Q_noCSIT);
                                
                                
                                
                                
                                
                                
                                

                                
                                
                                
                                                             
                                




                end
                
                
                Capacity_CSIT(i) = Temp ./ samples;
                Capacity_noCSIT(i) = Temp1 ./ samples;
                %Simulation ends here...  
                
                
                [SNR_dB'      Capacity_CSIT       Capacity_noCSIT]



end


% Plot of Waterfilling capacity (simulation) versus SNR
plot(SNR_dB,Capacity_CSIT);
hold on;
plot(SNR_dB,Capacity_noCSIT);
hold on;



    

    
    