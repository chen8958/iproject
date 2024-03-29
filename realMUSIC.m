function [angel_plant]=realMUSIC(MicPos)
c=343;
path=cd;

%%
% 
% 
[D MicNum]=size(MicPos);
    %% Windowing
% cd([path '\audio_R'])
% [x1 fs]=audioread('female_16k_10s.wav');
% SorLen=fs*4;
% Source=[x1(1:fs*4)];
% cd(path)
fs=16000;
%%
% [s1 fs]=audioread('p1.wav');
% [s2 fs]=audioread('p2.wav');
% s=[fft(s1).';fft(s2).'];
%%
% load('p1.mat');
% s1=m;
% load('p2.mat');
% s2=m;
% s=[fft(s1);fft(s2)];
%%
w=wgn(1000,1,0);
load('p1.mat');
%s1=conv(m,w);
s1=m;
s1=fft(s1);
for i=1:length(s1)
if abs(s1(i))>1
    s1(i)=s1(i)/abs(s1(i));
end
if abs(s1(i))<0.5
    s1(i)=s1(i)/abs(s1(i));
end
end
% s1=m;
%s1=s1/a;
load('p2.mat');
%s2=conv(m,w);
s2=m;

s2=fft(s2);
for i=1:length(s2)
if abs(s2(i))>1
    s2(i)=s2(i)/abs(s2(i));
end
if abs(s2(i))<0.5
    s2(i)=s2(i)/abs(s2(i));
end
end
% s2=m;
%s2=s2/a;
s=[s1;s2];

%%

NFFT=length(s);
df=fs/NFFT;
Freqs=0:df:(NFFT/2-1)*df;
Rxx=zeros(MicNum,MicNum,length(Freqs));


for ff=1:1:length(Freqs)
    Rxx_tmp(:,:,ff) = s(:,ff)*s(:,ff)';
end
Rxx=Rxx+Rxx_tmp;

Ps=zeros(MicNum,MicNum,length(Freqs));
P=zeros(MicNum,MicNum,length(Freqs));
for ff=1:1:length(Freqs);
[V,D]=eig(Rxx(:,:,ff));
[D_sort,index] = sort(diag(D),'descend');
V_sort = V(:,index);
for SignalSub=1:2
Ps(:,:,ff)=Ps(:,:,ff)+(V_sort(:,SignalSub)*V_sort(:,SignalSub)');
end
P(:,:,ff)=eye(MicNum)-Ps(:,:,ff);
end



    angel_curve=zeros(1,360);
    angel_plant=zeros(length(Freqs),360);
    for deg=1:360
        for ff=1:length(Freqs)

            k = 2*pi*Freqs(ff)/c; 
            kappa=[cosd(deg),sind(deg),0];
            for MicNo=1:MicNum
                a(MicNo,1)=exp(1i*k*kappa*MicPos(:,MicNo));
            end
%             w=a';
%             w=(inv(Rxx(:,:,ff)+0.01*eye(MicNum))*a)/(a'*inv(Rxx(:,:,ff)+0.01*eye(MicNum))*a);
%             angel_curve(deg)=angel_curve(deg)+abs(w*s(:,ff));
%             angel_plant(ff,deg)=abs(w*s(:,ff));
            angel_curve(deg)=angel_curve(deg)+abs(1/(a'*((P(:,:,ff)+0.00001*eye(2))*a)));
            angel_plant(ff,deg)=abs(1/(a'*((P(:,:,ff)+0.00001*eye(2))*a)));
        end
    end
    
    figure(1)
    deg=1:360;
    %pcolor(deg,Freqs,angel_plant);
    contourf(deg,Freqs,angel_plant);
    shading interp;
    figure(2)
    plot(deg,angel_curve);
    
end