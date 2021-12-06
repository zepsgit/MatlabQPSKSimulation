clear all
close all
M = 4;
span=10;
sps=10;
data_number = 4e4; 

data1 = randi([0 M-1],data_number,1);
data2 = randi([0 M-1],data_number,1);
data3 = randi([0 M-1],data_number,1);
data4 = randi([0 M-1],data_number,1);

txSig1 = pskmod(data1,M,pi/M);
txSig2 = pskmod(data2,M,pi/M);
txSig3 = pskmod(data3,M,pi/M);
txSig4 = pskmod(data4,M,pi/M);

txSigL1=txSig1;
txSigL2=[txSig1,txSig2]/sqrt(2);
txSigL4=[txSig1,txSig2,txSig3,txSig4]/sqrt(4);
txSig20db_n=awgn(txSigL1,20,'measured');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% investigate the constellation with different SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for snr_db=-10:5:15
        txSig_n = awgn(txSigL1,snr_db,'measured');
        scatterplot(txSig_n);
        title(['SNR ' num2str(snr_db) ' dB']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%investigate the psd with different rolloff factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for rolloff_factor=0.1:0.3:1
    Root_raised_cosine_filter = rcosdesign(rolloff_factor, span, sps, 'sqrt');
    upconverted_Transmitted_signal = upfirdn(txSig20db_n, Root_raised_cosine_filter, sps);
    [pxx frequency]=pwelch(upconverted_Transmitted_signal,hamming(1024),[],[],data_number*sps/2,'centered');
    plot(frequency,10*log10(pxx))
    legend('Alpha = 0.10','Alpha = 0.40','Alpha = 0.70','Alpha = 1.00');
    ylabel('PSD (dB/Hz)') 
    xlabel('Frequency (Hz)') 
    hold on
    grid on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%investigate the implication of an error in phase estimate
% Also can use phasenoise = comm.PhaseNoise('Level',[-70 -80])
%       phasenoise(txSig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_error_base1=pi/6;
phase_error_base2=pi/8;
r = 2*rand(1,1)-1;% r in range of [-1 1]
phase_errors1=r*phase_error_base1;% random phase in [-pi/4,pi/4]
phase_errors2=r*phase_error_base2;% random phase in [-pi/8,pi/8]
% plot the the constellation with random phase error
figure;
txSig1_phase_err1 = pskmod(data1,M,phase_errors1);
txSig1_phase_err2 = pskmod(data1,M,phase_errors2);
EbNo = (0:0.5:10)';% awgn should use symbol energy No ratio, need to covert
ber_unest1=zeros(length(EbNo),1);
ber_est1=zeros(length(EbNo),1);

ber_unest2=zeros(length(EbNo),1);
ber_est2=zeros(length(EbNo),1);
berQ = berawgn(EbNo,'psk',M,'diff');

for i=1:length(EbNo)
    txSig1_phase_err_n_1 = awgn(txSig1_phase_err1,10*log10(2)+EbNo(i),'measured');
    txSig1_phase_err_n_2 = awgn(txSig1_phase_err2,10*log10(2)+EbNo(i),'measured');
    phase_est_1=angle(mean(txSig1_phase_err_n_1.^4))/4;
    phase_est_2=angle(mean(txSig1_phase_err_n_2.^4))/4;
    demod_txSig1_phase_err_n_1=pskdemod(txSig1_phase_err_n_1, M, 0);
    demod_txSig1_phase_err_n_est_1=pskdemod(txSig1_phase_err_n_1,M,phase_est_1);

    demod_txSig1_phase_err_n_2=pskdemod(txSig1_phase_err_n_2, M, 0);
    demod_txSig1_phase_err_n_est_2=pskdemod(txSig1_phase_err_n_2,M,phase_est_2);
   
    numErrs_unest_1 = symerr(data1,demod_txSig1_phase_err_n_1);
    numErrs_est_1=symerr(data1,demod_txSig1_phase_err_n_est_1);
    ber_unest1(i)=numErrs_unest_1/data_number;
    ber_est1(i)=numErrs_est_1/data_number;
%{
For phase error in [-pi/8 pi/8]
    numErrs_unest_2 = symerr(data1,demod_txSig1_phase_err_n_1);
    numErrs_est_2=symerr(data1,demod_txSig1_phase_err_n_est_1);
    ber_unest2(i)=numErrs_unest_2/data_number;
    ber_est2(i)=numErrs_est_2/data_number;
%}
end
h1=semilogy(EbNo,berQ, EbNo,ber_unest1,'ro', EbNo,ber_est1,'go');
set(h1(1),'linewidth',1);
set(h1(2),'linewidth',1);
set(h1(3),'linewidth',1);

xlabel('Eb/No (dB)')
ylabel('BER')
legend('AWGN theory', 'PhaseErr1 in [-pi/6 pi/6]','Phase est for PhaseErr1');
title("Theoretical BER vs Real");
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% caculate bit error probability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
EbNo = (0:0.5:10)';% awgn should use symbol energy No ratio, need to covert
ber1=zeros(length(EbNo),1);
ber1Ray=zeros(length(EbNo),1);
ber2=zeros(length(EbNo),1);
ber2Ray=zeros(length(EbNo),1);
ber4=zeros(length(EbNo),1);
ber4Ray=zeros(length(EbNo),1);
berQ = berawgn(EbNo,'psk',M,'diff');
RayberQ = berfading(EbNo,'psk',M,1);
%SNR=10*log10(2)+EbNo;
%RayberQ=0.5*(1-1./sqrt(1+1./SNR));
R=raylrnd(1,data_number,1);
RSigL1=txSigL1.*R;
RSigL2=txSigL2.*R;
RSigL4=txSigL4.*R;
for i=1:length(EbNo)
    txSigL1_n = awgn(txSigL1,10*log10(2)+EbNo(i),'measured');
    % construct rayleigh channel using two independent gaussian channel
    %railey noise consist of two equally gaussion, so noise divided into 2
    %for gaussian channel, s/n=2Eb/No; for railey s/n=2Eb/(No/2)
    demod_txSigL1_n=pskdemod(txSigL1_n, M, pi/4);
    numErrs1 = symerr(data1,demod_txSigL1_n);
    ber1(i)=numErrs1/data_number;

    raySigL1=awgn(RSigL1,10*log10(4)+EbNo(i),'measured');    
    demod_raySigL1=pskdemod(raySigL1, M, pi/4);
    numErrs1Ray = symerr(data1,demod_raySigL1);
    ber1Ray(i)=numErrs1Ray/data_number;

    txSigL2_n = awgn(txSigL2,10*log10(2)+EbNo(i));
    demod_txSigL2_n=pskdemod(txSigL2_n, M, pi/4);
    raySigL2=awgn(RSigL2,10*log10(4)+EbNo(i),'measured');
    demod_raySigL2=pskdemod(raySigL2, M, pi/4);
    numErrs2Ray_1 = symerr(data1,demod_raySigL2(:,1));
    numErrs2Ray_2 = symerr(data2,demod_raySigL2(:,2));
    numErrs2_1 = symerr(data1,demod_txSigL2_n(:,1));
    numErrs2_2 = symerr(data2,demod_txSigL2_n(:,2));
    numErrs2=numErrs2_1+numErrs2_2;
    ber2(i)=numErrs2/(data_number*2);
    numErrs2Ray = numErrs2Ray_1+numErrs2Ray_2;
    ber2Ray(i)=numErrs2Ray/(2*data_number);

    txSigL4_n = awgn(txSigL4,10*log10(2)+EbNo(i));
    demod_txSigL4_n=pskdemod(txSigL4_n, M, pi/4);
    numErrs4_1 = symerr(data1,demod_txSigL4_n(:,1));
    numErrs4_2 = symerr(data2,demod_txSigL4_n(:,2));
    numErrs4_3 = symerr(data3,demod_txSigL4_n(:,3));
    numErrs4_4 = symerr(data4,demod_txSigL4_n(:,4));
    numErrs4=numErrs4_1+numErrs4_2+numErrs4_3+numErrs4_4;
    ber4(i)=numErrs4/(data_number*4);

    raySigL4=awgn(RSigL4,10*log10(4)+EbNo(i),'measured');    
    demod_raySigL4=pskdemod(raySigL4, M, pi/4);
    numErrs4Ray_1 = symerr(data1,demod_raySigL4(:,1));
    numErrs4Ray_2 = symerr(data2,demod_raySigL4(:,2));
    numErrs4Ray_3 = symerr(data3,demod_raySigL4(:,3));
    numErrs4Ray_4 = symerr(data4,demod_raySigL4(:,4));
    numErrs4Ray=numErrs4Ray_1+numErrs4Ray_2+numErrs4Ray_3+numErrs4Ray_4;
    ber4Ray(i)=numErrs4Ray/(4*data_number);
end
h=semilogy(EbNo,berQ,EbNo, ber1,'g-.o',EbNo, ber2,'k:x',EbNo, ber4,'m-s',EbNo,RayberQ,'r',EbNo,ber1Ray,'g*',EbNo,ber2Ray,'kx',EbNo,ber4Ray,'ms');
set(h(1),'linewidth',3);
set(h(2),'linewidth',1);
set(h(3),'linewidth',1);
set(h(4),'linewidth',1);
set(h(5),'linewidth',3);
set(h(6),'linewidth',1);
set(h(7),'linewidth',1);
set(h(8),'linewidth',1);
xlabel('Eb/No (dB)')
ylabel('BER')
legend('AWGN theory','L=1','L=2', 'L=4','Ray theory','Ray real L=1','Ray real L=2','Ray real L=4');
title("Theoretical BER vs Real");
grid on;
