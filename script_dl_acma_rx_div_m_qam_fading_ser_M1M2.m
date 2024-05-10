%% %%%%%%%%%%%%%%        Adaptive Constellation Multiple Access %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Dr Indu L Shakya, 22 Feb 2023 %%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [simSer, Thpt, simSer_, Thpt_,opt_rot] = script_dl_acma_rx_div_m_qam_fading_ser_M1M2(esno_1,esno_2,n,m1,m2,p1,p2,l,ch,rotate)

% [simSer, Thpt, simSer_, Thpt_]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2(20,20,1e+4,16,16,0.9,0.1,1,1,0)
%clear
rp=[];
N = n; % number of symbols
M1 = m1; % number of constellation points source 1 
M2 = m2;  % number of constellation points source 2
T = 2;

L = l;
nErr=0;
nErr_=0;

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

% Optimal phase rotation for User 1 for minimum Euclidian distance
Vlen=20;
Rot_Vec=zeros(1,Vlen);
SNR_loss_M12_full=zeros(1,Vlen);
sq_dist_col=zeros(1,Vlen);
sq_dist_M12_full=zeros(1,Vlen);

if rotate==100
    100;
    for q=1: length(Rot_Vec)
    [SNR_loss_M12_full(q),sq_dist_col(q),sq_dist_M12_full(q)]=CSMA_minimum_dist_pa_analysis(M1,M2,p1,p2,(2*pi)/(Vlen-q+1),0);
    Rot_Vec(q)=2*pi/(Vlen-q+1);
    end


[max_d, ind_rot]=max(sq_dist_col);

rotat=Rot_Vec(ind_rot);

opt_rot=rotat;

end




%%Rotated Conestellation User 1
if rotate
MQAM1_const_mag=abs(MQAM1_const);
MQAM1_const_ang=angle(MQAM1_const);
rotat=rotate;
MQAM1_const=MQAM1_const_mag.*exp(i*(MQAM1_const_ang+rotat));
%scatterplot(MQAM1_const)
end


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
    
    coll_MQAM=[MQAM1_const.' repmat(MQAM2_const(1,b),length(MQAM1_const),1)];
    colla_MQAM= [ colla_MQAM; coll_MQAM];
end

%scatterplot(colla_MQAM(:,1)+colla_MQAM(:,2))


%ipHat = zeros(1,N); % init
%ipHat_ = zeros(1,N); % init

%for ii = 1:5:length(Es_N0_dB)

    P1=p1;
    P2=p2;

for n_sym=1:N
    
    ip_orig = randsrc(1,1,alphaMqam1) + j*randsrc(1,1,alphaMqam1);
    ip=ip_orig;
    
    if rotate
        ip_const_mag=abs(ip);
        ip_const_ang=angle(ip);
        
        rotat1=rotat;
        ip=ip_const_mag.*exp(i*(ip_const_ang+rotat1));
        
    end
    
    ip_= randsrc(1,1,alphaMqam2) + j*randsrc(1,1,alphaMqam2);
    
    %scatterplot(ip)
    
    d_poss=zeros(M1*M2,1);
    
    s = repmat(k*ip,L,1); % normalization of energy to 1
    s_= repmat(k_*ip_,L,1);
    
    
    n = 1/sqrt(2)*[randn(L,1) + j*randn(L,1)]; % white guassian noise, 0dB variance
    n__ = 1/sqrt(2)*[randn(L,1) + j*randn(L,1)]; % white guassian noise, 0dB variance
    
    
    if ch==0
    h= 1/sqrt(2)*[randn(L,1) + j*randn(L,1)];    
    h__= 1/sqrt(2)*[randn(L,1) + j*randn(L,1)];
    
    else
        
    h = ones(L,1);
    h__= ones (L,1);
    end
    

    r = h.*(s*sqrt(P1) + s_*sqrt(P2))+ 10^(-Es_N0_1_dB/20)*n; % additive white gaussian noise
    r__ = h__.*(s*sqrt(P1) + s_*sqrt(P2))+ 10^(-Es_N0_2_dB/20)*n__; % additive white gaussian noise
    
    rp=[rp r];
    
    
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
   
 %   colla_MQAM(ind,1)
    
    ipHat=reshape(colla_MQAM(ind,1),1,length(ip));
    
%      if rotate
%         ipHat_const_mag=abs(ipHat);
%         ipHat_const_ang=angle(ipHat);
%         
%         rotat1n=rotate;
%         ipHat=ipHat_const_mag.*exp(i*(ipHat_const_ang+rotat1n));
%      end
    
    %ip_Hat=reshape(colla_MQAM(ind,2),1,length(ip_));
    
        nErr=nErr+size(find([ip- ipHat]),2);
    
    %nErr = size(find([ip- ipHat]),2); % counting the number of errors
    %nErr_ = size(find([ip_- ip_Hat]),2); % counting the number of errors
    

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
    nErr_ = nErr_+size(find([ip_- ip_Hat_]),2); % counting the number of errors
    
end    


simSer = nErr/N;
Thpt= log2(M1)*(1-simSer);


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