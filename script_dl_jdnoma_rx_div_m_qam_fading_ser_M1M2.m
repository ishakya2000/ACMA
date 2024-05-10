%% Collaborative Joint Detection for LTE-A heretegenous networks

function [simSer, Thpt, simSer_, Thpt_] = script_dl_colla_rx_div_m_qam_fading_ser_M1M2n(esno_1,esno_2,n,m1,m2,p1,p2,l,ch)

% [simSer, Thpt, simSer_, Thpt_]= script_dl_colla_rx_div_m_qam_fading_ser_M1M2n(20,20,1e+4,16,16,0.9,0.1,2,0)
%clear
N = n; % number of symbols
M1 = m1; % number of constellation points source 1 
M2 = m2;  % number of constellation points source 2
T = 2;

L = l;

%ratio_lin=10^(ratio/10);

if M1==2 
    MQAM1_const=(2*de2bi(0:(2^M1-1))-1)*(-1);
end

if M2==2 
    MQAM2_const=(2*de2bi(0:(2^M2-1))-1)*(-1);
end


k = sqrt(1/((2/3)*(M1-1))); % normalizing factor
k_ = sqrt(1/((2/3)*(M2-1))); % normalizing factor


m1 = [1:sqrt(M1)/2]; % alphabets
alphaMqam1 = [-(2*m1-1) 2*m1-1]; 

m2 = [1:sqrt(M2)/2]; % alphabets
alphaMqam2 = [-(2*m2-1) 2*m2-1]; 

Es_N0_1_dB = esno_1; %[0:40]; % multiple Es/N0 values
Es_N0_2_dB = esno_2; %[0:40]; % multiple Es/N0 values


for n=1:length(alphaMqam1)
    for m=1:length(alphaMqam1)
        MQAM1(n,m)=alphaMqam1(1,m)+j*alphaMqam1(1,n);
    end
end

MQAM1_const=reshape(MQAM1,1,m*n);

clear n,m;

for n=1:length(alphaMqam2)
    for m=1:length(alphaMqam2)
        MQAM2(n,m)=alphaMqam2(1,m)+j*alphaMqam2(1,n);
    end
end

MQAM2_const=reshape(MQAM2,1,m*n);

clear n,m;

colla_MQAM=[];
for b=1:length(MQAM2_const)
    
    coll_MQAM=[MQAM1_const' repmat(MQAM2_const(1,b),length(MQAM1_const),1)];
    colla_MQAM= [ colla_MQAM; coll_MQAM];
end

%scatterplot(MQAM_const)


ipHat = zeros(1,N); % init
ipHat_ = zeros(1,N); % init

%for ii = 1:5:length(Es_N0_dB)

    ip = randsrc(1,N,alphaMqam1) + j*randsrc(1,N,alphaMqam1);
    ip_= randsrc(1,N,alphaMqam2) + j*randsrc(1,N,alphaMqam2);
    
    d_poss=zeros(M1*M2,N);
    
    s = repmat(k*ip,L,1); % normalization of energy to 1
    s_= repmat(k_*ip_,L,1);
    
    P1=p1;
    P2=p2;
    
    n = 1/sqrt(2)*[randn(L,N) + j*randn(L,N)]; % white guassian noise, 0dB variance
    n__ = 1/sqrt(2)*[randn(L,N) + j*randn(L,N)]; % white guassian noise, 0dB variance
    
    
    if ch==0
    h= 1/sqrt(2)*[randn(L,N) + j*randn(L,N)];    
    h__= 1/sqrt(2)*[randn(L,N) + j*randn(L,N)];
    
    else
        
    h = ones(L,N);
    h__= ones (L,N);
    end
    

    r = h.*(s*sqrt(P1) + s_*sqrt(P2))+ 10^(-Es_N0_1_dB/20)*n; % additive white gaussian noise
    r__ = h__.*(s*sqrt(P1) + s_*sqrt(P2))+ 10^(-Es_N0_2_dB/20)*n__; % additive white gaussian noise
    
    
    r_=[];
    for n=1:L
    r_=[r_; repmat(r(n,:),M1*M2,1)];
    end
    
    clear r;
        

    if L==1        
        
        d_poss = abs(r_(1:M1*M2,:)-(kron(h(1,:),(sqrt(P1)*k*colla_MQAM(:,1)+sqrt(P2)*k_*colla_MQAM(:,2))))).^2;
         
    else
        
    d_poss = abs(r_(1:M1*M2,:)-(kron(h(1,:),(sqrt(P1)*k*colla_MQAM(:,1)+sqrt(P2)*k_*colla_MQAM(:,2))))).^2 ...
             + abs(r_(M1*M2+1:M1*M2*L,:)-(kron(h(2,:),(sqrt(P1)*k*colla_MQAM(:,1)+sqrt(P2)*k_*colla_MQAM(:,2))))).^2;
    end
    
    [dmin,ind]=min(d_poss);
    
    ipHat=reshape(colla_MQAM(ind,1),1,length(ip));
    %ip_Hat=reshape(colla_MQAM(ind,2),1,length(ip_));
    
    
    nErr = size(find([ip- ipHat]),2); % counting the number of errors
    %nErr_ = size(find([ip_- ip_Hat]),2); % counting the number of errors
    



simSer = nErr/N;
Thpt= log2(M1)*(1-simSer);








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%Detection of weak %%%%%%%%%%%%%%%%%%%%signal%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



r___=[];
    for n=1:L
    r___=[r___; repmat(r__(n,:),M1*M2,1)];
    end
    
    clear r__;
        

    if L==1        
        
        d_poss = abs(r___(1:M1*M2,:)-(kron(h__(1,:),(sqrt(P1)*k*colla_MQAM(:,1)+sqrt(P2)*k_*colla_MQAM(:,2))))).^2;
         
    else
        
    d_poss = abs(r___(1:M1*M2,:)-(kron(h__(1,:),(sqrt(P1)*k*colla_MQAM(:,1)+sqrt(P2)*k_*colla_MQAM(:,2))))).^2 ...
             + abs(r___(M1*M2+1:M1*M2*L,:)-(kron(h__(2,:),(sqrt(P1)*k*colla_MQAM(:,1)+sqrt(P2)*k_*colla_MQAM(:,2))))).^2;
    end
    
    [dmin,ind]=min(d_poss);
    
    %ipHat=reshape(colla_MQAM(ind,1),1,length(ip));
    ip_Hat_=reshape(colla_MQAM(ind,2),1,length(ip_));
    
    
    %nErr = size(find([ip- ipHat]),2); % counting the number of errors
    nErr_ = size(find([ip_- ip_Hat_]),2); % counting the number of errors










simSer_ = nErr_/N;
Thpt_= log2(M2)*(1-simSer_);
%theorySer = 2*(1-1/sqrt(M))*erfc(k*sqrt((10.^(Es_N0_dB/10)))) ...
 %             - (1-2/sqrt(M) + 1/M)*(erfc(k*sqrt((10.^(Es_N0_dB/10))))).^2;
%close all
%figure
%semilogy(Es_N0_dB,theorySer,'bs-','LineWidth',2);
%hold on
%semilogy(Es_N0_dB,simSer,'m*-','Linewidth',1);
%hold on;
%axis([0 40 10^-5 1])
%grid on
%legend('theory', 'simulation');
%xlabel('Es/No, dB')
%ylabel('Symbol Error Rate')
%title('Symbol error probability curve for M-QAM modulation')





