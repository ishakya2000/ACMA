function [SNR_loss_M12_full,sq_dist_col,sq_dist_M12_full]=CSMA_minimum_dist_pa_analysis(M1,M2,p1,p2,rotate,sh_const)

%[SNR_loss_M12_full,sq_dist_col,sq_dist_M12_full]=CSMA_minimum_dist_pa_analysis(16,4,0.8333,0.1677,0)

%M1 = 16; % number of constellation points source 1 
%M2 = 4;  % number of constellation points source 2
%rotate=1;

if M1==2 
    MQAM1_const=(2*de2bi(0:(2^M1-1))-1)*(-1);
end

if M2==2 
    MQAM2_const=(2*de2bi(0:(2^M2-1))-1)*(-1);
end


k = sqrt(1/((2/3)*(M1-1))); % normalizing factor
k_ = sqrt(1/((2/3)*(M2-1))); % normalizing factor
kk_ = sqrt(1/((2/3)*(M1*M2-1))); % normalizing factor



m1 = [1:sqrt(M1)/2]; % alphabets
alphaMqam1 = [-(2*m1-1) 2*m1-1]; 

m2 = [1:sqrt(M2)/2]; % alphabets
alphaMqam2 = [-(2*m2-1) 2*m2-1]; 

m12=[1:sqrt(M1*M2)/2]; % alphabets
alphaMqam12 = [-(2*m12-1) 2*m12-1]; 


for n=1:length(alphaMqam1)
    for m=1:length(alphaMqam1)
        MQAM1(n,m)=alphaMqam1(1,m)+j*alphaMqam1(1,n);
    end
end

MQAM1_const=reshape(MQAM1,1,m*n);

%%Rotated Conestellation
if rotate
MQAM1_const_mag=abs(MQAM1_const);
MQAM1_const_ang=angle(MQAM1_const);
rotat=rotate;
MQAM1_const=MQAM1_const_mag.*exp(i*(MQAM1_const_ang+rotat));
%scatterplot(MQAM1_const_rotated)
end


clear n,m;

for n=1:length(alphaMqam2)
    for m=1:length(alphaMqam2)
        MQAM2(n,m)=alphaMqam2(1,m)+j*alphaMqam2(1,n);
    end
end

MQAM2_const=reshape(MQAM2,1,m*n);

clear n,m;


for n=1:length(alphaMqam12)
    for m=1:length(alphaMqam12)
        MQAM12(n,m)=alphaMqam12(1,m)+j*alphaMqam12(1,n);
    end
end

MQAM12_const=reshape(MQAM12,1,m*n);

clear n,m;

colla_MQAM=[];
for b=1:length(MQAM2_const)
    
    coll_MQAM=[MQAM1_const' repmat(MQAM2_const(1,b),length(MQAM1_const),1)];
    colla_MQAM= [ colla_MQAM; coll_MQAM];
end

%scatterplot(MQAM1_const)




for n=1:size(MQAM1_const,2)
for m=1:size(MQAM2_const,2)
joint(m,n)=sqrt(p1)*k*MQAM1_const(n)+sqrt(p2)*k_*MQAM2_const(m);
end
end
colla_const=reshape(joint,M1*M2,1);

M1_const=sqrt(p1)*k*MQAM1_const;
%scatterplot(M1_const)

M2_const=sqrt(p2)*k_*MQAM2_const;
%scatterplot(M2_const)

M12_const=sqrt(p1+p2)*kk_*MQAM12_const;

if sh_const

    scatterplot(colla_const);
    scatterplot(M12_const);

end


for e=1:length(M1_const)
    for f=1:length(M1_const)
        dist_M1(e,f)=abs(M1_const(e)-M1_const(f)).^2;
    end
end

sq_dist_M1 = min(dist_M1(dist_M1 > 0));

clear e,f;

for e=1:length(M12_const)
    for f=1:length(M12_const)
        dist_M12(e,f)=abs(M12_const(e)-M12_const(f)).^2;
    end
end


idx = eye(length(M12_const),length(M12_const));
%sq_dist_M12 = min(dist_M12(dist_M12 > 0))

dist_M12_= (1-idx).*dist_M12;
dist_M12__=dist_M12_(~idx);
sq_dist_M12 = min(dist_M12__);


for b=1:length(colla_const)
    for c=1:length(colla_const)
        dist_colla(b,c)=abs(colla_const(b)-colla_const(c)).^2;
    end
end

dist_colla_=(1-idx).*dist_colla;
dist_colla__=dist_colla_(~idx);
%sq_dist_col= min(dist_colla(dist_colla > 0));
sq_dist_col= min(min(dist_colla__));
avg_dist_col=mean(mean(dist_colla__));



M1_const_full=k*MQAM1_const;
%scatterplot(M1_const_full)

M12_const_full=kk_*MQAM12_const;
%scatterplot(M1_const_full)

for e=1:length(M1_const_full)
    for f=1:length(M1_const_full)
        dist_M1_full(e,f)=abs(M1_const_full(e)-M1_const_full(f)).^2;
    end
end


%sq_dist_M1_full = min(dist_M1_full(dist_M1_full > 0));
sq_dist_M1_full = min(min(dist_M1_full));


clear e,f;


for e=1:length(M12_const_full)
    for f=1:length(M12_const_full)
        dist_M12_full(e,f)=abs(M12_const_full(e)-M12_const_full(f)).^2;
    end
end

dist_M12_full_=(1-idx).*dist_M12_full;
dist_M12_full__=dist_M12_full_(~idx);
%sq_dist_M12_full = min(dist_M12_full(dist_M12_full > 0));
sq_dist_M12_full = min(min(dist_M12_full__));



%SNR_loss_M1= 10*log10(sq_dist_M1/sq_dist_col)
%SNR_loss_M1_full= 10*log10(sq_dist_M1_full/sq_dist_col)

if (sq_dist_M12==0) | (sq_dist_col==0)
    SNR_loss_M12=100;
    SNR_loss_M12_full=100;
else
     
SNR_loss_M12= 10*log10(sq_dist_M12/sq_dist_col);
SNR_loss_M12_full= 10*log10(sq_dist_M12_full/sq_dist_col);

end

graph=0;

if graph
    
ps=p1;
pw=p2;
simSer=zeros(1,16);

for n=1:16
    snr(n)=(n-1)*3;
    %[simSer(n)] = compute_symbol_error_rate(snr(n), M1*M2); %AWGN channels
    simSer(n) = script_m_qam_fading_ser(snr(n),100000,M1*M2); %Fading channels
end

figure;
semilogy(snr,simSer,'ro');

hold on

offsets=10*log10(ps/1);
offsetw=10*log10(pw/1);

semilogy(snr+(SNR_loss_M12_full),simSer,'bo');


for n=1:16
    snr(n)=(n-1)*3;
    %[simSer(n)] = compute_symbol_error_rate(snr(n), M1*M2); %AWGN channels
    simSerm2(n) = script_m_qam_fading_ser(snr(n),100000,M1); %Fading channels
end

hold on

offsets=10*log10(ps/1);
offsetw=10*log10(pw/1);

%[sermqamfading, j_sermqamfading] = joint_detection_ser_analysis(M1,M2,snr-offsets)
%semilogy(snr-(offsetw),sermqamfading, 'g*');


% Lower bound on SER of weak user U2
semilogy(snr-(offsetw),simSerm2,'b*');


 hold on
% Upper bound on SER of strong user U1

semilogy(snr-(offsets),simSerm2,'b>');


for n=1:16
    snr(n)=(n-1)*3;
    %[simSer(n)] = compute_symbol_error_rate(snr(n), M1*M2); %AWGN channels
    %simSerm2(n) = script_m_qam_fading_ser(snr(n),100000,M1); %Fading channels
    ana_ser(n)=mqam_ser_rayleigh(snr(n),M1*M2);

end

hold on


semilogy(snr+(SNR_loss_M12_full),ana_ser,'go');
end

