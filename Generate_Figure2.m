%% Generate Figure 2 of I. L. Shakya and F. H. Ali, "Adaptive Constellation Multiple Access for Beyond 5G Wireless Systems," 
%% in IEEE Wireless Communications Letters, doi: 10.1109/LWC.2024.3367924


clear;
rotation =[]
rot=[];
for n=1:20
[simSer_sa4, Thpt_sa4, simSer_wa4, Thpt_wa4,rota]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2(30,30,1e+3,4,4,n/20,1-(n/20),1,0,100);
rot=[rot rota];
end
rotation = [rotation;rot];

figure;
plot(1/20:1/20:1,(rot*360)./(2*pi),'-+');
%legend()
hold on

rot=[];
for n=1:20
[simSer_sa4, Thpt_sa4, simSer_wa4, Thpt_wa4,rota]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2(30,30,1e+3,4,16,n/20,1-(n/20),1,0,100);
rot=[rot rota];
end
rotation = [rotation;rot];

plot(1/20:1/20:1,(rot*360)./(2*pi),'r-+');
%legend('ACMA 4-QAM/16-QAM')
hold on;


rot=[];
for n=1:20
[simSer_sa4, Thpt_sa4, simSer_wa4, Thpt_wa4,rota]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2(30,30,1e+3,16,4,n/20,1-(n/20),1,0,100);
rot=[rot rota];
end
rotation = [rotation;rot];


plot(1/20:1/20:1,(rot*360)./(2*pi),'g-+');
%legend('ACMA 16-QAM/4-QAM');
hold on;


rot=[];
for n=1:20
[simSer_sa4, Thpt_sa4, simSer_wa4, Thpt_wa4,rota]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2(30,30,1e+3,16,16,n/20,1-(n/20),1,0,100);
rot=[rot rota];
end
rotation = [rotation;rot];

hold on;
grid on;
grid minor;

plot(1/20:1/20:1,(rot*360)./(2*pi),'y-+');
legend('ACMA 4-QAM/4-QAM', 'ACMA 4-QAM/16-QAM','ACMA 16-QAM/4-QAM', 'ACMA 16-QAM/16-QAM');

ylim([0 200]);
xlabel('Power allocation factor \alpha_{1}, \alpha_{2}=1-\alpha_{1}');
ylabel('Phase rotation of U1 \delta_{1}, degrees');


%rotdeg=rotation*360/(2*pi)
%pow=1/20:1/20:1
%%%%

%grid on
