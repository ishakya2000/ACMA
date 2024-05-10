
%% Generate Figure 5 of I. L. Shakya and F. H. Ali, "Adaptive Constellation Multiple Access for Beyond 5G Wireless Systems," 
%% in IEEE Wireless Communications Letters, May 2024, doi: 10.1109/LWC.2024.3367924


%16-4QAM
clear;
%pack;
p2=0.5:0.05:0.95;
p1=1-p2;

Ser_s_sic=zeros(1,length(p2));
Ser_w_sic=zeros(1,length(p2));
Ser_s_col=zeros(1,length(p2));
Ser_w_col=zeros(1,length(p2));

Ser_s_ac=zeros(1,length(p2));
Ser_w_ac=zeros(1,length(p2));

Run_num=1e+5;

l=100;
for a=1:l
    a
for n=1:length(p2)
    
[simSer_s_sic(n), Thpt_s_sic(n), simSer_w_sic(n), Thpt_w_sic(n)] =script_dl_sic_rx_div_m_qam_fading_ser_M1M2(40,40,Run_num,16,4,p1(n),p2(n),1,0);

[simSer_s_col(n), Thpt_s_col(n), simSer_w_col(n), Thpt_w_col(n)]= script_dl_jdnoma_rx_div_m_qam_fading_ser_M1M2(40,40,Run_num,16,4,p1(n),p2(n),1,0);

[simSer_sa(n), Thpt_sa(n), simSer_wa(n), Thpt_wa(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2(40,40,Run_num,16,4,p1(n),p2(n),1,0,100);

end
Ser_s_sic=Ser_s_sic+simSer_s_sic;
Ser_w_sic=Ser_w_sic+simSer_w_sic;
Ser_s_col=Ser_s_col+simSer_s_col;
Ser_w_col=Ser_w_col+simSer_w_col;
Ser_s_ac=Ser_s_ac+simSer_sa;
Ser_w_ac=Ser_w_ac+simSer_wa;
end

simSer_s_sic_m=Ser_s_sic/l;
simSer_w_sic_m=Ser_w_sic/l;
simSer_s_col_m=Ser_s_col/l;
simSer_w_col_m=Ser_w_col/l;

simSer_s_sa_m=Ser_s_ac/l;
simSer_w_sa_m=Ser_w_ac/l;


figure;
semilogy(p1,simSer_s_sic_m,'b-o')
hold on;
semilogy(p1,simSer_w_sic_m,'b-+')
hold on;
semilogy(p1,simSer_s_col_m,'r-o')
hold on;
semilogy(p1,simSer_w_col_m,'r-+')
hold on;
semilogy(p1,simSer_s_sa_m,'g-o')
hold on;
semilogy(p1,simSer_w_sa_m,'g-+')

legend('PD-NOMA U1','PD-NOMA U2','JD-NOMA U1','JD-NOMA U2','ACMA U1','ACMA U2');


%%%%16-16QAM

clear;
%pack;
p2=0.5:0.05:0.95;
p1=1-p2;

Ser_s_sic=zeros(1,length(p2));
Ser_w_sic=zeros(1,length(p2));
Ser_s_col=zeros(1,length(p2));
Ser_w_col=zeros(1,length(p2));

Ser_s_ac=zeros(1,length(p2));
Ser_w_ac=zeros(1,length(p2));

Run_num=1e+5;

l=100;
for a=1:l
    a
for n=1:length(p2)
    
[simSer_s_sic(n), Thpt_s_sic(n), simSer_w_sic(n), Thpt_w_sic(n)] =script_dl_sic_rx_div_m_qam_fading_ser_M1M2(50,50,Run_num,16,16,p1(n),p2(n),1,0);

[simSer_s_col(n), Thpt_s_col(n), simSer_w_col(n), Thpt_w_col(n)]= script_dl_jdnoma_rx_div_m_qam_fading_ser_M1M2(50,50,Run_num,16,16,p1(n),p2(n),1,0);

[simSer_sa(n), Thpt_sa(n), simSer_wa(n), Thpt_wa(n)]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2(50,50,Run_num,16,16,p1(n),p2(n),1,0,100);

end
Ser_s_sic=Ser_s_sic+simSer_s_sic;
Ser_w_sic=Ser_w_sic+simSer_w_sic;
Ser_s_col=Ser_s_col+simSer_s_col;
Ser_w_col=Ser_w_col+simSer_w_col;
Ser_s_ac=Ser_s_ac+simSer_sa;
Ser_w_ac=Ser_w_ac+simSer_wa;
end

simSer_s_sic_m=Ser_s_sic/l;
simSer_w_sic_m=Ser_w_sic/l;
simSer_s_col_m=Ser_s_col/l;
simSer_w_col_m=Ser_w_col/l;

simSer_s_sa_m=Ser_s_ac/l;
simSer_w_sa_m=Ser_w_ac/l;


figure;
semilogy(p1,simSer_s_sic_m,'b-o')
hold on;
semilogy(p1,simSer_w_sic_m,'b-+')
hold on;
semilogy(p1,simSer_s_col_m,'r-o')
hold on;
semilogy(p1,simSer_w_col_m,'r-+')
hold on;
semilogy(p1,simSer_s_sa_m,'g-o')
hold on;
semilogy(p1,simSer_w_sa_m,'g-+')

legend('PD-NOMA U1','PD-NOMA U2','JD-NOMA U1','JD-NOMA U2','ACMA U1','ACMA U2')