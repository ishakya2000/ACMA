%% Superposition Transmssion with SIC

function [simSer_s, Thpt_s, simSer_f, Thpt_f]= script_dl_sic_rx_div_m_qam_fading_ser_M1M2(esno_1,esno_2,n,m1,m2,p1,p2,l,ch)

% [simSer_s, Thpt_s, simSer_f, Thpt_f] =script_dl_sic_rx_div_m_qam_fading_ser_M1M2(20,0,1e+6,4,4,0.5,0.5,1,0)
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
    n_ = 1/sqrt(2)*[randn(L,N) + j*randn(L,N)]; % white guassian noise, 0dB variance
    
    
    if ch==0
    h= 1/sqrt(2)*[randn(L,N) + j*randn(L,N)];    
    h_= 1/sqrt(2)*[randn(L,N) + j*randn(L,N)];
    
    else
        
    h = ones(L,N);
    h_= ones (L,N);
    end
    

   
    r = h.*(s*sqrt(P1) + s_*sqrt(P2))+ 10^(-Es_N0_1_dB/20)*n; % additive white gaussian noise      
    r_ = h_.*(s*sqrt(P1) + s_*sqrt(P2))+ 10^(-Es_N0_2_dB/20)*n_; % additive white gaussian noise
    
    
    %r_=[];
    %for n=1:L
    %r_=[r_; repmat(r(n,:),M1*M2,1)];
    %end
    
    %clear r;
        

    %%%%%%%%%%%%%%%%%%%%%%%Detection of strong users%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     % demodulation
    %y_=r_./h_;
    yw=sum(conj(h).*r,1)./sum(h.*conj(h),1);
    
    yw_re = real(yw)/(k_*sqrt(P2)); % real part
    yw_im = imag(yw)/(k_*sqrt(P2)); % imaginary part
    
    
    
    % rounding to the nearest alphabet
    % 0 to 2 --> 1
    % 2 to 4 --> 3
    % 4 to 6 --> 5 etc
    
    ipHat_wre = 2*floor(yw_re/2)+1;
    ipHat_wre(find(ipHat_wre>max(alphaMqam2))) = max(alphaMqam2);
    ipHat_wre(find(ipHat_wre<min(alphaMqam2))) = min(alphaMqam2);
	     
    % rounding to the nearest alphabet
    % 0 to 2 --> 1
    % 2 to 4 --> 3
    % 4 to 6 --> 5 etc
    ipHat_wim = 2*floor(yw_im/2)+1;
    ipHat_wim(find(ipHat_wim>max(alphaMqam2))) = max(alphaMqam2);
    ipHat_wim(find(ipHat_wim<min(alphaMqam2))) = min(alphaMqam2);
    
    
    ip_Hat_w = ipHat_wre + j*ipHat_wim;
    
    
    
    nErr_w = size(find([ip_- ip_Hat_w]),2); 
    
    simSer_w = nErr_w/N;
    
    
    
    r_s=r-(h.*(sqrt(P2)*k_*repmat(ip_Hat_w,L,1)));
    %y__=r__./h_;
    y_s=sum(conj(h).*r_s,1)./sum(h.*conj(h),1);
    
    y_sre = real(y_s)/(k*sqrt(P1)); % real part
    y_sim = imag(y_s)/(k*sqrt(P1)); % imaginary part
    
    
    
    % rounding to the nearest alphabet
    % 0 to 2 --> 1
    % 2 to 4 --> 3
    % 4 to 6 --> 5 etc
    
    ipHat_sre = 2*floor(y_sre/2)+1;
    ipHat_sre(find(ipHat_sre>max(alphaMqam1))) = max(alphaMqam1);
    ipHat_sre(find(ipHat_sre<min(alphaMqam1))) = min(alphaMqam1);
	     
    % rounding to the nearest alphabet
    % 0 to 2 --> 1
    % 2 to 4 --> 3
    % 4 to 6 --> 5 etc
    ipHat_sim = 2*floor(y_sim/2)+1;
    ipHat_sim(find(ipHat_sim>max(alphaMqam1))) = max(alphaMqam1);
    ipHat_sim(find(ipHat_sim<min(alphaMqam1))) = min(alphaMqam1);
    
    
    ipHat_s = ipHat_sre + j*ipHat_sim; 
    
    nErrs = size(find([ip- ipHat_s]),2); 
    
    
    simSer_s = nErrs/N;
    
    Thpt_s= log2(M1)*(1-simSer_s);
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Detection of weak%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%signal%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    %y=r./h;
    y_=sum(conj(h_).*r_,1)./sum(h_.*conj(h_),1);
    
    y_re = real(y_)/(k_*sqrt(P2)); % real part
    y_im = imag(y_)/(k_*sqrt(P2)); % imaginary part
    
    
    
    % rounding to the nearest alphabet
    % 0 to 2 --> 1
    % 2 to 4 --> 3
    % 4 to 6 --> 5 etc
    
    ipHat_re = 2*floor(y_re/2)+1;
    ipHat_re(find(ipHat_re>max(alphaMqam2))) = max(alphaMqam2);
    ipHat_re(find(ipHat_re<min(alphaMqam2))) = min(alphaMqam2);
	     
    % rounding to the nearest alphabet
    % 0 to 2 --> 1
    % 2 to 4 --> 3
    % 4 to 6 --> 5 etc
    ipHat_im = 2*floor(y_im/2)+1;
    ipHat_im(find(ipHat_im>max(alphaMqam2))) = max(alphaMqam2);
    ipHat_im(find(ipHat_im<min(alphaMqam2))) = min(alphaMqam2);
    
    
    ipHat = ipHat_re + j*ipHat_im; 
    
    nErr_f = size(find([ip_- ipHat]),2); 
    
    
    simSer_f = nErr_f/N;
    
    Thpt_f= log2(M2)*(1-simSer_f);
    
  
    
