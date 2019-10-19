cd(['.\audio_R'])
[x1 fs]=audioread('female_16k_10s.wav');
[x2 fs]=audioread('male_16k_10s.wav');
SorLen=fs*4;
Source=[x1(1:fs*4);x2(1:fs*4)];
cd('..');


load('freq_azi75mic0.mat');
m=resample(m,16000,44100);
sor45mic0=m;
load('freq_azi75mic15.mat');
m=resample(m,16000,44100);
sor45mic1=m;
load('freq_azi285mic0.mat');
m=resample(m,16000,44100);
sor135mic0=m;
load('freq_azi285mic15.mat');
m=resample(m,16000,44100);
sor135mic1=m;
%%
%45 female x1 135 male x2
sor45mic0=conv(sor45mic0,x1);
sor45mic1=conv(sor45mic1,x1);
sor135mic0=conv(sor135mic0,x2);
sor135mic1=conv(sor135mic1,x2);

m=sor45mic0+sor135mic0;
m=m';
save(['mp1.mat'],'m');

m=sor45mic1+sor135mic1;
m=m';
save(['mp2.mat'],'m');