%% Generate Figure 4 of I. L. Shakya and F. H. Ali, "Adaptive Constellation Multiple Access for Beyond 5G Wireless Systems," 
%% in IEEE Wireless Communications Letters, May 2024, doi: 10.1109/LWC.2024.3367924%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Figure 6
%%
%%% Throughput from SER of TDMA against ACMA 16-16 QAM, K=2

%PD-NOMA-SIC%%%%%%%%%%%%%%
clear;
N=40;
NumRun=1e+5;

%%%%%%%%%%%%%%%%%%%%%%TDMA equal power and equal channel gains

% N=40;
% NumRun=10e+2;
P1=0.5
P2=0.5
M1=16
M2=16
Off=0

% PD-NOMA
for n=1:N
[simSer_s_sic1(n), Thpt_s_sic1(n), simSer_w_sic1(n), Thpt_w_sic1(n)] =script_dl_sic_rx_div_m_qam_fading_ser_M1M2((n-1)*2,((n-1)*2-Off),NumRun,M1,M2,P1,P2,1,0);
end
%figure
hold on
plot(0:2:(N-1)*2, Thpt_s_sic1+Thpt_w_sic1,'g-^');


%JD-NOMA%%%%%%%%%%%%%%

for n=1:N
[simSer_sc1(n), Thpt_sc1(n), simSer_wc1(n), Thpt_wc1(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2((n-1)*2,((n-1)*2-Off),NumRun,M1,M2,P1,P2,1,0,0);
end

hold on
plot(0:2:(N-1)*2, Thpt_sc1+Thpt_wc1,'b-^');

%ACMA%%%%%%%%%%%%%%

%N=40
for n=1:N
[simSer_sa1(n), Thpt_sa1(n), simSer_wa1(n), Thpt_wa1(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2((n-1)*2,((n-1)*2-Off),NumRun,M1,M2,P1,P2,1,0,100);
end

hold on
plot(0:2:(N-1)*2, Thpt_sa1+Thpt_wa1,'r-+');

%TDMA
for n=1:N
[simSer_tdma_s_u1(n), Thpt_tdma_s_u1(n)] = script_m_qam_fading_ser(((n-1)*2),NumRun,M1*M2);
[simSer_tdma_w_u2(n), Thpt_tdma_w_u2(n)] = script_m_qam_fading_ser(((n-1)*2)-Off,NumRun,M1*M2);
[simSer_single(n), Thpt_single(n)] = script_m_qam_fading_ser(((n-1)*2),NumRun,M1*M2);
end

Sum_Thput_tdma=(P1*Thpt_tdma_s_u1)+(P2*Thpt_tdma_w_u2);
%figure;
hold on;
plot(0:2:(N-1)*2, Sum_Thput_tdma,'y-o');

%% Unequal power and 20 dB SNr offset

P1=0.1
P2=0.9
M1=16
M2=16
Off=20
for n=1:N
[simSer_s_sic1(n), Thpt_s_sic1(n), simSer_w_sic1(n), Thpt_w_sic1(n)] =script_dl_sic_rx_div_m_qam_fading_ser_M1M2((n-1)*2,((n-1)*2-Off),NumRun,M1,M2,P1,P2,1,0);
end

%figure
hold on
plot(0:2:(N-1)*2, Thpt_s_sic1+Thpt_w_sic1,'g^');

%JD-NOMA%%%%%%%%%%%%%%

for n=1:N
[simSer_sc1(n), Thpt_sc1(n), simSer_wc1(n), Thpt_wc1(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2((n-1)*2,((n-1)*2-Off),NumRun,M1,M2,P1,P2,1,0,0);
end

hold on
plot(0:2:(N-1)*2, Thpt_sc1+Thpt_wc1,'b^');

%ACMA%%%%%%%%%%%%%%

%N=40
for n=1:N
[simSer_sa1(n), Thpt_sa1(n), simSer_wa1(n), Thpt_wa1(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2((n-1)*2,((n-1)*2-Off),NumRun,M1,M2,P1,P2,1,0,100);
end

hold on
plot(0:2:(N-1)*2, Thpt_sa1+Thpt_wa1,'r+');


% P1=0.1; 
% P2=0.9;
% N=40
for n=1:N
[simSer_tdma_s_u1(n), Thpt_tdma_s_u1(n)] = script_m_qam_fading_ser(((n-1)*2),NumRun,M1*M2);
[simSer_tdma_w_u2(n), Thpt_tdma_w_u2(n)] = script_m_qam_fading_ser(((n-1)*2-Off),NumRun,M1*M2);
[simSer_single(n), Thpt_single(n)] = script_m_qam_fading_ser(((n-1)*2),NumRun,M1*M2);
end

Sum_Thput_tdma=(P1*Thpt_tdma_s_u1)+(P2*Thpt_tdma_w_u2);
%figure;
hold on;
plot(0:2:(N-1)*2, Sum_Thput_tdma,'y+');


for n=1:N
[simSer_sing_256(n), Thpt_sing_256(n)] = script_m_qam_fading_ser((n-1)*2,1e+5,256)

end

hold on;
plot(0:2:(N-1)*2, Thpt_sing_256,'mo');


legend('PD-NOMA 0.5,0.5','JD-NOMA 0.5,0.5','ACMA 0.5,0.5','TDMA 0.5,0.5','PD-NOMA 0.1,0.9','JD-NOMA 0.1,0.9','ACMA 0.1,0.9','TDMA 0.1,0.9', 'Single-user 256-QAM')
