The Matlab codes for Figure 2-6 used in the paper I. L. Shakya and F. H. Ali, "Adaptive Constellation Multiple Access for Beyond 5G Wireless Systems," 
%% in IEEE Wireless Communications Letters, May 2024, doi: 10.1109/LWC.2024.3367924

Notes:
- Codes have been tested for Matlab 2019a
- After downloading the files, please store all in one folder. Run Matlab and change the path to that folder. Simply type: Generate Figure2, Generate Figure3, and so on. Codes may take few hours (Figures 2,3,4,6) or two days to complete for Figure 5 (you can see progress). 
- The codes may run on earlier versions of Matlab or Octave too. But there is no guarantee that they will work. 
  
- Figure 1 of the paper can be generate by operning the file script_dl_acma_rx_div_m_qam_fading_ser_M1M2.m and removing comment % from Line 106 and saving the file.
  - %scatterplot(colla_MQAM(:,1)+colla_MQAM(:,2))
  - After than please run the line below. Say we want to see composite constellation of JD-NOMA that rotates the second user's constellation 9.4 degrees,  
      [simSer, Thpt, simSer_, Thpt_]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2(20,20,1e+4,16,16,0.9,0.1,1,1,pi/(180/9.4))
  -For ACMA it automatically find the best constellation. Just enter 100 for the rotation value to show the constellation:
      [simSer, Thpt, simSer_, Thpt_]= script_dl_acma_rx_div_m_qam_fading_ser_M1M2(20,20,1e+4,16,16,0.9,0.1,1,1,100)
   - Add % in the above file to stop it from disrupting code runs (you will not need too many scatterplots popping out) for Figures 3-6

- If you fiind this work (paper and codes) useful, please acknowledge by citing the work. Thank you. Indu L. Shakya
