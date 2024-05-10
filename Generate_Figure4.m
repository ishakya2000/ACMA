%% Generate Figure 4 of I. L. Shakya and F. H. Ali, "Adaptive Constellation Multiple Access for Beyond 5G Wireless Systems," 
%% in IEEE Wireless Communications Letters, May 2024, doi: 10.1109/LWC.2024.3367924%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Figure 4 a)

clear;
Run_num=1e+5;

% Non-equal power
P1=0.75;

for n=1:11
[simSer_s_a(n), Thpt_s_a(n), simSer_w_a(n), Thpt_w_a(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,16,P1,1-P1,1,0,100);
%[simSer_s_sic(n), Thpt_s_sic(n), simSer_w_sic(n), Thpt_w_sic(n)] =script_dl_sic_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,16,P1,1-P1,1,0);
[simSer_s_col(n), Thpt_s_col(n), simSer_w_col(n), Thpt_w_col(n)]= script_dl_jdnoma_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,16,P1,1-P1,1,0);
[simSer_s_ad(n), Thpt_s_ad(n), simSer_w_ad(n), Thpt_w_ad(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,16,P1,1-P1,2,0,100);
end

figure;
semilogy(20:5:50,simSer_s_a(5:11),'g-o');
hold on;
semilogy(20:5:50,simSer_w_a(5:11),'g-+');
hold on;
semilogy(20:5:50,simSer_s_col(5:11),'r-o');
hold on;
semilogy(20:5:50,simSer_w_col(5:11),'r-+');
hold on;
semilogy(20:5:50,simSer_s_ad(5:11),'b-o');
hold on;
semilogy(20:5:50,simSer_w_ad(5:11),'b-+');

for n=1:25
[simSer_sing_256(n), Thpt_sing_256(n)] = script_m_qam_fading_ser((n-1)*2,1e+5,256)

end

hold on
semilogy(20:2:48, simSer_sing_256(11:25),'y^');


legend('ACMA U1','ACMA U2','JD-NOMA U1','JD-NOMA U2','ACMA U1 Dual','ACMA U2 Dual','Single User 256-QAM')


%% Figure 4 b) 

% Equal power
P1=0.5;
Run_num=1e+4;


for n=1:11
[simSer_s_a(n), Thpt_s_a(n), simSer_w_a(n), Thpt_w_a(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,16,P1,1-P1,1,0,100);
[simSer_s_sic(n), Thpt_s_sic(n), simSer_w_sic(n), Thpt_w_sic(n)] =script_dl_sic_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,16,P1,1-P1,1,0);
[simSer_s_col(n), Thpt_s_col(n), simSer_w_col(n), Thpt_w_col(n)]= script_dl_jdnoma_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,16,P1,1-P1,1,0);
[simSer_s_ad(n), Thpt_s_ad(n), simSer_w_ad(n), Thpt_w_ad(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,16,P1,1-P1,2,0,100);
end

figure;
semilogy(20:5:50,simSer_s_a(5:11),'g-o');
hold on;
semilogy(20:5:50,simSer_w_a(5:11),'g-+');
hold on;
semilogy(20:5:50,simSer_s_col(5:11),'r-o');
hold on;
semilogy(20:5:50,simSer_w_col(5:11),'r-+');
hold on;
semilogy(20:5:50,simSer_s_sic(5:11),'m-o');
hold on;
semilogy(20:5:50,simSer_w_sic(5:11),'m-+');
hold on;
semilogy(20:5:50,simSer_s_ad(5:11),'b-o');
hold on;
semilogy(20:5:50,simSer_w_ad(5:11),'b-+');

for n=1:25
[simSer_sing_256(n), Thpt_sing_256(n)] = script_m_qam_fading_ser((n-1)*2,1e+5,256)

end

hold on
semilogy(20:2:48, simSer_sing_256(11:25),'y^');

legend('ACMA U1','ACMA U2','JD-NOMA U1','JD-NOMA U2','PD-NOMA U1','PD-NOMA U2','ACMA U1 Dual','ACMA U2 Dual','Single User 256-QAM')
