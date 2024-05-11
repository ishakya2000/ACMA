

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [simSer, Thpt] = script_m_qam_fading_ser(esno,n,m1)

% Creative Commons
% Attribution-Noncommercial 2.5 India
% You are free:
% to Share — to copy, distribute and transmit the work
% to Remix — to adapt the work
% Under the following conditions:
% Attribution. You must attribute the work in the manner 
% specified by the author or licensor (but not in any way 
% that suggests that they endorse you or your use of the work). 
% Noncommercial. You may not use this work for commercial purposes. 
% For any reuse or distribution, you must make clear to others the 
% license terms of this work. The best way to do this is with a 
% link to this web page.
% Any of the above conditions can be waived if you get permission 
% from the copyright holder.
% Nothing in this license impairs or restricts the author's moral rights.
% http://creativecommons.org/licenses/by-nc/2.5/in/

% Script for simulating M-QAM (4-QAM, 16-QAM, 64-QAM, 256QAM, 1024QAM) 
% transmission and reception and compare the simulated and theoretical 
% symbol error probability

% Checked for proper operation with Octave Version 3.0.0
% Author	: Krishna Pillai
% Email		: krishna@dsplog.com
% Version	: 1.0
% Date		: 28 April 2008
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear
N =  n; % number of symbols
M = m1; % number of constellation points      

k = sqrt(1/((2/3)*(M-1))); % normalizing factor

m = [1:sqrt(M)/2]; % alphabets
alphaMqam = [-(2*m-1) 2*m-1]; 

Es_N0_dB = esno; %[0:5:40]; % multiple Es/N0 values


for n=1:length(alphaMqam)
    for m=1:length(alphaMqam)
        MQAM(n,m)=alphaMqam(1,m)+j*alphaMqam(1,n);
    end
end

MQAM_const=reshape(MQAM,1,m*n);
%scatterplot(MQAM_const)

nRx =  1;

ipHat = zeros(1,N); % init

jj = nRx;

    %for ii = 1:length(Es_N0_dB)

        ip = randsrc(1,N,alphaMqam) + j*randsrc(1,N,alphaMqam);
        s = k*ip; % normalization of energy to 1
        
        n = 1/sqrt(2)*[randn(nRx,N) + j*randn(nRx,N)]; % white gaussian noise, 0dB variance
        h = 1/sqrt(2)*[randn(nRx,N) + j*randn(nRx,N)]; % Rayleigh channel

        sD = kron(ones(nRx,1),s);
        y = h.*sD + 10^(-Es_N0_dB/20)*n;

        
        % equalization maximal ratio combining 
        yHat =  sum(conj(h).*y,1)./sum(h.*conj(h),1); 

        y_re = real(yHat)/k; % real part
        y_im = imag(yHat)/k; % imaginary part



        % rounding to the nearest alphabet
        % 0 to 2 --> 1
        % 2 to 4 --> 3
        % 4 to 6 --> 5 etc
        ipHat_re = 2*floor(y_re/2)+1;
        ipHat_re(find(ipHat_re>max(alphaMqam))) = max(alphaMqam);
        ipHat_re(find(ipHat_re<min(alphaMqam))) = min(alphaMqam);

        % rounding to the nearest alphabet
        % 0 to 2 --> 1
        % 2 to 4 --> 3
        % 4 to 6 --> 5 etc
        ipHat_im = 2*floor(y_im/2)+1;
        ipHat_im(find(ipHat_im>max(alphaMqam))) = max(alphaMqam);
        ipHat_im(find(ipHat_im<min(alphaMqam))) = min(alphaMqam);

        ipHat = ipHat_re + j*ipHat_im; 
        nErr = size(find([ip- ipHat]),2); % counting the number of errors
    

    %end
    
simSer = nErr/N;
Thpt= log2(M)*(1-simSer);

end

% %theorySer = 2*(1-1/sqrt(M))*erfc(k*sqrt((10.^(Es_N0_dB/10)))) ...
%  %             - (1-2/sqrt(M) + 1/M)*(erfc(k*sqrt((10.^(Es_N0_dB/10))))).^2;
% %close all
% %figure
% %semilogy(Es_N0_dB,theorySer,'bs-','LineWidth',2);
% hold on
% semilogy(Es_N0_dB,simSer,'m*-','Linewidth',1);
% hold on;
% %axis([0 40 10^-5 1])
% grid on
% %legend('theory', 'simulation');
% xlabel('Es/No, dB')
% ylabel('Symbol Error Rate')
% title('Symbol error probability curve for M-QAM modulation with nRx antenna')
% 




