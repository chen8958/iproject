function [angel_plant]=realDAS(MicPos)
c=343;
path=cd;

%%
% 
% 
[D MicNum]=size(MicPos);
    %% Windowing
cd([path '\audio_R'])
[x1 fs]=audioread('female_16k_10s.wav');
SorLen=fs*4;
Source=[x1(1:fs*4)];
cd(path)
% % %%
% [s1 fs]=audioread('p1.wav');
% [s2 fs]=audioread('p2.wav');
% s=[fft(s1).';fft(s2).'];
%%
%white noise

% w=wgn(1000,1,0);
% load('p1.mat');
% %s1=conv(m,w);
% s1=m;
% %s1=s1/a;
% load('p2.mat');
% %s2=conv(m,w);
% s2=m;
% %s2=s2/a;
% s=[fft(s1);fft(s2)];

%%

% 
load('rp1.mat');
s1=m;
load('rp2.mat');
s2=m;
s=[fft(s1);fft(s2)];

%%
%windowing
% NWIN=1024;
% hopsize=NWIN/2;                                                            % 50% overlap
% NumOfFrame=2*floor(fs*4/NWIN)-1;                                           % number of frames
% win = hann(NWIN+1);                                                        % hanning window
% win = win(1:end-1).';
% %% FFT
% NFFT=2^nextpow2(NWIN);
% df=fs/NFFT;
% Freqs=0:df:(NFFT/2-1)*df;
NFFT=length(s);
df=fs/NFFT;
Freqs=0:df:(NFFT/2-1)*df;


    angel_curve=zeros(1,360);
    angel_plant=zeros(length(Freqs),360);
    for deg=1:360
        for ff=1:length(Freqs)

            k = 2*pi*Freqs(ff)/c; 
            kappa=[cosd(deg),sind(deg),0];
            for MicNo=1:MicNum
                a(MicNo,1)=exp(1i*k*kappa*MicPos(:,MicNo));
            end
            w=a';
            angel_curve(deg)=angel_curve(deg)+abs(w*s(:,ff));
            angel_plant(ff,deg)=abs(w*s(:,ff));
        end
    end
   
    
    figure(1);
    deg=1:360;
    pcolor(deg,Freqs,angel_plant);
%     contourf(deg,Freqs,angel_plant);
    shading interp;
    figure(2);
    plot(deg,angel_curve);
    %%
   %dip
     for deg=1:360
        for ff=1:length(Freqs)
            angel_plant(ff,deg)=angel_plant(ff,deg)^2;
        end
     end
    figure(3);
    deg=1:360;
    pcolor(deg,Freqs,angel_plant);
%     contourf(deg,Freqs,angel_plant);
    shading interp;
    
end