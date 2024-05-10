
%% Generate Figure 3 of I. L. Shakya and F. H. Ali, "Adaptive Constellation Multiple Access for Beyond 5G Wireless Systems," 
%% in IEEE Wireless Communications Letters, May 2024, doi: 10.1109/LWC.2024.3367924

%Figure 3 a)

clear;
Run_num=1e+5;

% Non-equal power
P1=0.1;

for n=1:11
%[simSer(b,n), Thpt(b,n), simSer_(b,n), Thpt_(b,n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2((n-1)*2,(n-1)*2,1e+4,M_((b-1)*2+1),M_((b-1)*2+2),0.5,0.5,1,1,100);
[simSer_s_a(n), Thpt_s_a(n), simSer_w_a(n), Thpt_w_a(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,4,P1,1-P1,1,0,100);
[simSer_s_sic(n), Thpt_s_sic(n), simSer_w_sic(n), Thpt_w_sic(n)] =script_dl_sic_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,4,P1,1-P1,1,0);
[simSer_s_col(n), Thpt_s_col(n), simSer_w_col(n), Thpt_w_col(n)]= script_dl_jdnoma_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,4,P1,1-P1,1,0);
[simSer(n), Thpt(n), simSer_(n), Thpt_(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,4,P1,1-P1,1,0,pi/(180/9.4));
end

figure;
semilogy(0:5:50,simSer_s_a,'g-o');
hold on;
semilogy(0:5:50,simSer_w_a,'g-+');
hold on;
semilogy(0:5:50,simSer_s_col,'r-o');
hold on;
semilogy(0:5:50,simSer_w_col,'r-+');
hold on;
semilogy(0:5:50,simSer_s_sic,'b-o');
hold on;
semilogy(0:5:50,simSer_w_sic,'b-+');
hold on;
semilogy(0:5:50,simSer,'-o');
hold on;
semilogy(0:5:50,simSer_,'-+');

legend('ACMA U1','ACMA U2','JD-NOMA U1 0^{o}','JD-NOMA U2 0^{o}','PD-NOMA U1','PD-NOMA U2','JD-NOMA U1 9.4^{o}','JD-NOMA U2 9.4^{o}')


%% Figure 3 b) with different phase offsets

P1=0.35;

for n=1:11
[simSer_s_a(n), Thpt_s_a(n), simSer_w_a(n), Thpt_w_a(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,4,P1,1-P1,1,0,100);
[simSer_s_sic(n), Thpt_s_sic(n), simSer_w_sic(n), Thpt_w_sic(n)] =script_dl_sic_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,4,P1,1-P1,1,0);
[simSer_s_col(n), Thpt_s_col(n), simSer_w_col(n), Thpt_w_col(n)]= script_dl_jdnoma_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,4,P1,1-P1,1,0);
[simSer_s_chan(n), Thpt_s_chan(n), simSer_w_chan(n), Thpt_w_chan(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2_mml((n-1)*5,(n-1)*5,Run_num,16,4,P1,1-P1,1,0,pi/(180/9.4),1);
[simSer(n), Thpt(n), simSer_(n), Thpt_(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2((n-1)*5,(n-1)*5,Run_num,16,4,P1,1-P1,1,0,pi/(180/9.4));
end
figure;
semilogy(0:5:50,simSer_s_a,'g-o');
hold on;
semilogy(0:5:50,simSer_w_a,'g-+');
hold on;
semilogy(0:5:50,simSer_s_col,'r-o');
hold on;
semilogy(0:5:50,simSer_w_col,'r-+');
hold on;
semilogy(0:5:50,simSer_s_sic,'b-o');
hold on;
semilogy(0:5:50,simSer_w_sic,'b-+');
hold on;
semilogy(0:5:50,simSer_s_chan,'y-o');
hold on;
semilogy(0:5:50,simSer_w_chan,'y-+');
hold on;
semilogy(0:5:50,simSer,'-o');
hold on;
semilogy(0:5:50,simSer_,'-+');

legend('ACMA U1','ACMA U2','JD-NOMA U1 0^{0}','JD-NOMA U2 0^{0}','PD-NOMA U1','PD-NOMA U2','JD-NOMA U1 9.4^{0}','JD-NOMA U2 9.4^{0}','JD-NOMA U1 2\pi/3','JD-NOMA U2 2\pi/3^{0}')
