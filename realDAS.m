function [angel_curve]=realDAS(MicPos)
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


% load('mp1.mat');
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
% load('mp2.mat');
% %s2=conv(m,w);
% s2=m;
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
% % 
% load('rp1.mat');
% s1=m;
% load('rp2.mat');
% s2=m;
% s=[fft(s1);fft(s2)];
%%
load('mp1.mat');
s1=m;
load('mp2.mat');
s2=m;
s=[fft(s1,1024);fft(s2,1024)];
% if(length(s1)>1024)
%     fot t
% end
NFFT=1024;
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
%             t=1:1/fs:length(Freqs);
%             sinusoid=fft(sin(2*pi*Freqs(ff)*t));
            w=a';
            angel_curve(deg)=angel_curve(deg)+abs(w*s(:,ff));
            angel_plant(ff,deg)=abs(w*s(:,ff));
%             angel_curve(deg)=angel_curve(deg)+abs(w*s(:,ff))*sinusoid(ff);
%             angel_plant(ff,deg)=abs(w*s(:,ff)*sinusoid(ff));
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
    
    %%
  %sep part
  [pks lo]=findpeaks(angel_curve);
  lo
  %lo=[15,165];
  s_all=[fft(s1);fft(s2)];
  
  NFFT=length(s_all);
  df=fs/NFFT;
  Freqs=0:df:(NFFT/2-1)*df;
  %%
  %set source number
  SorNum=2;
%   for ss=1:SorNum
%       kappa=[cosd(lo(ss)),sind(lo(ss)),0];
%       for ff=1:length(Freqs)
%           k = 2*pi*Freqs(ff)/c;
%           for MicNo=1:MicNum
%               a(MicNo,1)=exp(1i*k*kappa*MicPos(:,MicNo));
%           end
% %           w=a';
% %           w=inv(a'*a+0.001*eye(2))*a';
%           w=a'/(a'*a+0.001);
%           y_half(ss,ff)=w*s_all(:,ff);
%       end
%       y(ss,:)=[y_half(ss,:),zeros(1,1),fliplr(conj(y_half(ss,2:end)))];
%   end
      
    
      for ff=1:length(Freqs)
          k = 2*pi*Freqs(ff)/c;
          for ss=1:SorNum
              kappa=[cosd(lo(ss)),sind(lo(ss)),0];
          for MicNo=1:MicNum
              a(MicNo,ss)=exp(1i*k*kappa*MicPos(:,MicNo));
          end
          end
%           w=a';
%           w=inv(a'*a+0.001*eye(2))*a';
          w=inv(a'*a+0.001*eye(2))*a';
          y_half(:,ff)=w*s_all(:,ff);
      end
      for ss=1:SorNum
      y(ss,:)=[y_half(ss,:),zeros(1,1),fliplr(conj(y_half(ss,2:end)))];
      end
  
  
  for ss = 1:SorNum
      y(ss,:)=ifft(y(ss,:));
      audiowrite(['realTIKR_sep' num2str(MicNum) num2str(ss) '.wav'],y(ss,:)/max(abs(y(ss,:))),fs);
  end
    
end