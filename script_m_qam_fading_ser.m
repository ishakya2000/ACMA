

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [simSer, Thpt] = script_m_qam_fading_ser(esno,n,m1)

% Edited using codes by Krishna Pillai
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




