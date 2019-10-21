function [angel_curve]=realMVDR(MicPos)
c=343;
path=cd;

%%
% 
% 
[D MicNum]=size(MicPos);
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

% load('p1.mat');
% %s1=conv(m,w);
% s1=m;
% s1=fft(s1);
% for i=1:length(s1)
% if abs(s1(i))>1
%     s1(i)=s1(i)/abs(s1(i));
% end
% if abs(s1(i))<0.5
%     s1(i)=s1(i)/abs(s1(i));
% end
% end
% % s1=m;
% %s1=s1/a;
% load('p2.mat');
% %s2=conv(m,w);
% s2=m;
% 
% s2=fft(s2);
% for i=1:length(s2)
% if abs(s2(i))>1
%     s2(i)=s2(i)/abs(s2(i));
% end
% if abs(s2(i))<0.5
%     s2(i)=s2(i)/abs(s2(i));
% end
% end
% % s2=m;
% %s2=s2/a;
% s=[s1;s2];

%%
load('mp1.mat');
s1=m;
load('mp2.mat');
s2=m;
s=[fft(s1);fft(s2)];
%%

NFFT=length(s1);
df=fs/NFFT;
Freqs=0:df:(NFFT/2-1)*df;
Rxx=zeros(MicNum,MicNum,length(Freqs));


for ff=1:1:length(Freqs)
    Rxx_tmp(:,:,ff) = s(:,ff)*s(:,ff)';
end
Rxx=Rxx+Rxx_tmp;

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
            angel_curve(deg)=angel_curve(deg)+abs(1/(a'*inv(Rxx(:,:,ff)+0.01*eye(MicNum))*a));
            angel_plant(ff,deg)=abs(1/(a'*inv(Rxx(:,:,ff)+0.01*eye(MicNum))*a));
        end
    end
    
    figure(1)
    deg=1:360;
    %pcolor(deg,Freqs,angel_plant);
    contourf(deg,Freqs,angel_plant);
    shading interp;
    figure(2)
    plot(deg,angel_curve);
    
    %%
  %set source number
  [pks lo]=findpeaks(angel_curve);
  lo=[15,165];
  s_all=[fft(s1);fft(s2)];
  
  NFFT=length(s_all);
  df=fs/NFFT;
  Freqs=0:df:(NFFT/2-1)*df;
  
  SorNum=2;
  for ss=1:SorNum
      kappa=[cosd(lo(ss)),sind(lo(ss)),0];
      for ff=1:length(Freqs)
          k = 2*pi*Freqs(ff)/c;
          for MicNo=1:MicNum
              a(MicNo,1)=exp(1i*k*kappa*MicPos(:,MicNo));
          end
          w=inv(Rxx(:,:,ff)+0.01*eye(MicNum))*a/(a'*inv(Rxx(:,:,ff)+0.01*eye(MicNum))*a);
          y_half(ss,ff)=w'*s_all(:,ff);
      end
      y(ss,:)=[y_half(ss,:),zeros(1,1),fliplr(conj(y_half(ss,2:end)))];
  end
      
  for ss = 1:SorNum
      y(ss,:)=ifft(y(ss,:));
      audiowrite(['realMPDR_sep' num2str(MicNum) num2str(ss) '.wav'],y(ss,:)/max(abs(y(ss,:))),fs);
  end
    
end