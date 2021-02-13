%%% MATLAB code for AMC performance simulations
%%% -> an algorithm described in paper
%%% "Automatic Modulation Classification of Real Signals 
%%% in AWGN Channel Based on Sixth - Order Cumulants"
%%% M. Simic et al
%%%  - Radioengineering Journal, 2021

 
clc; clear all;
 
snr=10; % SNR value
N=250; % sample size N

% initiate numerical variables....
broj_2=0;
broj_4=0;
broj_1=0;
broj_8=0;
broj_16=0;
broj_32=0;
broj_64=0;
pogodak_c42=0;
pogodak_2_c42=0;
pogodak_4_c42=0;
pogodak_1_c42=0;
pogodak_8_c42=0;
pogodak_16_c42=0;
pogodak_32_c42=0;
pogodak_64_c42=0;
pogodak_c63=0;
pogodak_2_c63=0;
pogodak_4_c63=0;
pogodak_1_c63=0;
pogodak_8_c63=0;
pogodak_16_c63=0;
pogodak_32_c63=0;
pogodak_64_c63=0;
pogodak_c63_new=0;
pogodak_2_c63_new=0;
pogodak_4_c63_new=0;
pogodak_1_c63_new=0;
pogodak_8_c63_new=0;
pogodak_16_c63_new=0;
pogodak_32_c63_new=0;
pogodak_64_c63_new=0;

for brojac=1:2000
 
biraj=rand; % random modulation selection
if biraj>0.84
M = 4; % PAM-4
elseif biraj>0.70 M = 2; % BPSK
elseif biraj>0.56 M = 8; % PAM-8
elseif biraj>0.42 M = 16; % PAM-16
elseif biraj>0.28 M = 32; % PAM-32
elseif biraj>0.14 M = 64; % PAM-64
else M = 1; % QPSK
end;

if (M==1) % generate QPSK signal
    k = log2(4); 
    n = 4000*N/250; 
    nsamp = 1; 
 
    x = randint(n,1); 
    xsym = bi2de(reshape(x,k,length(x)/k).','left-msb');
    y1= modulate(modem.qammod(4),xsym);
    
else % generate real signal with modulation order M
    xr=randi([0 M-1],N,1);
    y1=pammod(xr,M);
end;
 
ytx = y1(1:N); % signal transmitted, sample size limit

h=1; % flat channel attenuation factor;
ytx=ytx*h;

snaga=0; % signal power estimation
for q=1:length(ytx)
    snaga=snaga+power(abs(ytx(q)),2);
end;
snaga=snaga/length(ytx);

ynoisy = awgn(ytx,snr,'measured'); % add noise
ssuma=snaga/power(10,snr/10); % noise power
 
yrx=ynoisy; % signal received

% normalized cumulant calculations:
cm4=0;
cm2=0;
c2=0;
cm6=0;
m41=0;
for q=1:length(yrx)
    cm4=cm4+power(abs(yrx(q)),4);
    c2=c2+power(yrx(q),2);
    cm2=cm2+power(abs(yrx(q)),2);
    cm6=cm6+power(abs(yrx(q)),6);    
    m41=m41+power((yrx(q)),2)*power(abs(yrx(q)),2);
end;
 
cm4=cm4/length(yrx);
cm2=cm2/length(yrx);
c2=c2/length(yrx);
cm6=cm6/length(yrx);
m41=m41/length(yrx);
 
c42=cm4-power(abs(c2),2)-2*power(cm2,2);
c42_n=c42/power((cm2-ssuma),2); % normalized 4th order cumulant value

c63=cm6-9*cm4*cm2+12*power(abs(c2),2)*cm2+12*power(cm2,3);
c63_n_old=c63/power((cm2-ssuma),3); % old normalized 4th order cumul. value

korekcija=6*power(abs(c2),2)*cm2-6*abs(c2)*abs(m41); 
korekcija=korekcija/power((cm2-ssuma),3); % normalized offset value

%%% decision making process - modulation classification
%%% 4th order cumulant:
  if c42_n<-1.68 odluka_c42=2; % decision: BPSK 
  elseif c42_n<-1.3 odluka_c42=4; % decision: PAM-4
  elseif c42_n<-1.224 odluka_c42=8; % decision: PAM-8
  elseif c42_n<-1.206 odluka_c42=16; % decision: PAM-16
  elseif c42_n<-1.201 odluka_c42=32; % decision: PAM-32
  elseif c42_n<-1.1 odluka_c42=64; % decision: PAM-64    
 else odluka_c42=1; % decision: QPSK
 end;
 if (odluka_c42==M) pogodak_c42=pogodak_c42+1; % count Pcc for C42...
 end;
 if ((odluka_c42==M)&&(M==2)) pogodak_2_c42=pogodak_2_c42+1;
 end;
 if ((odluka_c42==M)&&(M==4)) pogodak_4_c42=pogodak_4_c42+1;
 end;
 if ((odluka_c42==M)&&(M==1)) pogodak_1_c42=pogodak_1_c42+1;
 end;
 if ((odluka_c42==M)&&(M==8)) pogodak_8_c42=pogodak_8_c42+1;
 end;
 if ((odluka_c42==M)&&(M==16)) pogodak_16_c42=pogodak_16_c42+1;
 end;
 if ((odluka_c42==M)&&(M==32)) pogodak_32_c42=pogodak_32_c42+1;
 end;
 if ((odluka_c42==M)&&(M==64)) pogodak_64_c42=pogodak_64_c42+1;
 end;
 
%%% decision making process - modulation classification
%%% 6th order cumulant:
 if c63_n_old<7.83 odluka_c63=1; odluka_c63_new=1;% decision: QPSK
 
 else
     c63_n_new=c63_n_old+korekcija; % new 6th order cumulant - unbiased
      
     if  c63_n_new<6.87 odluka_c63_new=64; % decision: PAM-64
     elseif c63_n_new<6.91 odluka_c63_new=32; % decision: PAM-32
     elseif c63_n_new<7.06 odluka_c63_new=16; % decision: PAM-16    
     elseif c63_n_new<7.75 odluka_c63_new=8; % decision: PAM-8
     elseif c63_n_new<12.16 odluka_c63_new=4; % decision: PAM-4
     else odluka_c63_new=2; % decision: BPSK
     end;
 
 % decision with old 6th order cumulants...
     if  c63_n_old<11.66095 odluka_c63=64;
     elseif c63_n_old<11.67245 odluka_c63=32;
     elseif c63_n_old<11.72085 odluka_c63=16; 
     elseif c63_n_old<11.96 odluka_c63=8;
     elseif c63_n_old<14.08 odluka_c63=4;
     else odluka_c63=2;
     end;
 end;

 if (odluka_c63==M) pogodak_c63=pogodak_c63+1; % count Pcc for old C63...
 end;
 if ((odluka_c63==M)&&(M==2)) pogodak_2_c63=pogodak_2_c63+1;
 end;
 if ((odluka_c63==M)&&(M==4)) pogodak_4_c63=pogodak_4_c63+1;
 end;
 if ((odluka_c63==M)&&(M==1)) pogodak_1_c63=pogodak_1_c63+1;
 end;
 if ((odluka_c63==M)&&(M==8)) pogodak_8_c63=pogodak_8_c63+1;
 end;
 if ((odluka_c63==M)&&(M==16)) pogodak_16_c63=pogodak_16_c63+1;
 end;
 if ((odluka_c63==M)&&(M==32)) pogodak_32_c63=pogodak_32_c63+1;
 end;
 if ((odluka_c63==M)&&(M==64)) pogodak_64_c63=pogodak_64_c63+1;
 end;
 
 
 if (odluka_c63_new==M) pogodak_c63_new=pogodak_c63_new+1; % count Pcc for new C63...
 end;
 if ((odluka_c63_new==M)&&(M==2)) pogodak_2_c63_new=pogodak_2_c63_new+1;
 end;
 if ((odluka_c63_new==M)&&(M==4)) pogodak_4_c63_new=pogodak_4_c63_new+1;
 end;
 if ((odluka_c63_new==M)&&(M==1)) pogodak_1_c63_new=pogodak_1_c63_new+1;
 end;
 if ((odluka_c63_new==M)&&(M==8)) pogodak_8_c63_new=pogodak_8_c63_new+1;
 end;
 if ((odluka_c63_new==M)&&(M==16)) pogodak_16_c63_new=pogodak_16_c63_new+1;
 end;
 if ((odluka_c63_new==M)&&(M==32)) pogodak_32_c63_new=pogodak_32_c63_new+1;
 end;
 if ((odluka_c63_new==M)&&(M==64)) pogodak_64_c63_new=pogodak_64_c63_new+1;
 end;
 
 if (M==2) 
    broj_2=broj_2+1;
 elseif (M==4)
    broj_4=broj_4+1;
 elseif (M==8)
    broj_8=broj_8+1;
 elseif (M==16)
    broj_16=broj_16+1;
 elseif (M==32)
    broj_32=broj_32+1;
 elseif (M==64)
    broj_64=broj_64+1;
 else
    broj_1=broj_1+1; 
 end;
 
end;
 
pcc_c42_bpsk=pogodak_2_c42/broj_2; %Pcc of BPSK signal with C42
pcc_c42_4pam=pogodak_4_c42/broj_4; %Pcc of PAM-4 signal with C42
pcc_c42_8pam=pogodak_8_c42/broj_8; %Pcc of PAM-8 signal with C42
pcc_c42_16pam=pogodak_16_c42/broj_16; %Pcc of PAM-16 signal with C42
pcc_c42_32pam=pogodak_32_c42/broj_32; %Pcc of PAM-32 signal with C42
pcc_c42_64pam=pogodak_64_c42/broj_64; %Pcc of PAM-64 signal with C42
pcc_c42_qpsk=pogodak_1_c42/broj_1; %Pcc of QPSK signal with C42
PCC_c42=pogodak_c42/brojac %TOTAL PCC for 4th order cumulants
 
pcc_c63_bpsk_old=pogodak_2_c63/broj_2; %Pcc of BPSK signal with old C63
pcc_c63_4pam_old=pogodak_4_c63/broj_4; %Pcc of PAM-4 signal with old C63
pcc_c63_8pam_old=pogodak_8_c63/broj_8; %Pcc of PAM-8 signal with old C63
pcc_c63_16pam_old=pogodak_16_c63/broj_16; %Pcc of PAM-16 signal with old C63
pcc_c63_32pam_old=pogodak_32_c63/broj_32; %Pcc of PAM-32 signal with old C63
pcc_c63_64pam_old=pogodak_64_c63/broj_64; %Pcc of PAM-64 signal with old C63
pcc_c63_qpsk_old=pogodak_1_c63/broj_1; %Pcc of QPSK signal with old C63
PCC_c63_old=pogodak_c63/brojac %TOTAL PCC for old 6th order cumulants

pcc_c63_bpsk_new=pogodak_2_c63_new/broj_2; %Pcc of BPSK signal with new C63
pcc_c63_4pam_new=pogodak_4_c63_new/broj_4; %Pcc of PAM-4 signal with new C63
pcc_c63_8pam_new=pogodak_8_c63_new/broj_8; %Pcc of PAM-8 signal with new C63
pcc_c63_16pam_new=pogodak_16_c63_new/broj_16; %Pcc of PAM-16 signal with new C63
pcc_c63_32pam_new=pogodak_32_c63_new/broj_32; %Pcc of PAM-32 signal with new C63
pcc_c63_64pam_new=pogodak_64_c63_new/broj_64; %Pcc of PAM-64 signal with new C63
pcc_c63_qpsk_new=pogodak_1_c63_new/broj_1; %Pcc of QPSK signal with new C63
PCC_c63_new=pogodak_c63_new/brojac %TOTAL PCC for new 6th order cumulants
