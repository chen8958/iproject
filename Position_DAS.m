function [angel_plant]=Position_DAS(MicPos)
c=343;
SorNum=2;

[D MicNum]=size(MicPos);
    %% Windowing
      for i=1:MicNum
          [p(i,:) fs]=audioread("p"+i+".wav");
      end
%[p fs]=audioread("hello.wav");
%p=p';

NWIN=1024;
hopsize=NWIN/2;                                                            % 50% overlap
NumOfFrame=2*floor(fs*4/NWIN)-1;                                           % number of frames
win = hann(NWIN+1);                                                        % hanning window
win = win(1:end-1).';
%% FFT
NFFT=2^nextpow2(NWIN);
df=fs/NFFT;
Freqs=0:df:(NFFT/2-1)*df;

    for i=1:MicNum
    %x(i,1:4*fs)=p(i,1:4*fs);
    x(i,1:fs)=p(i,1:fs);
    fft_x(i,:)=fft(x(i,:),NFFT);
    end
    
    
    angel_curve=zeros(1,360);
    angel_plant=zeros(length(Freqs),360);
    for deg=1:360
        for ff=1:length(Freqs)
            k = 2*pi*Freqs(ff)/c; 
            kappa=[cosd(deg),sind(deg),0];
            for MicNo=1:MicNum
                a(MicNo)=exp(j*k*kappa*MicPos(:,MicNo));
            end
            w=conj(a);
            angel_curve(deg)=angel_curve(deg)+abs(w*fft_x(:,ff));
            angel_plant(ff,deg)=abs(w*fft_x(:,ff));
        end
    end
    figure(1)
    deg=1:360;
    pcolor(deg,Freqs,angel_plant);
    shading interp;
    figure(2)
    plot(deg,angel_curve);
    
    
    [argvalue,argmax] = max(angel_curve);

%     
% 
%     for FrameNo=1:NumOfFrame
%         % --time segment--
%         t_start=(FrameNo-1)*hopsize;
%         tt=(t_start+1):(t_start+NWIN);
% 
%         % --propagation--
%         for ff=1:length(Freqs)
%             k = 2*pi*Freqs(ff)/c;                                              
%             for ss = 1:SorNum
%             for m = 1:MicNum
%                  kappa=[cosd(argmax),sind(argmax),0];
%                  A(m,ss) =exp(1j*k*kappa*MicPos(:,m));
% %                  r = sqrt(sum((SorPos(:,ss)-MicPos(:,m)).^2));
% %                 A(m,ss) =exp(-1j*k*r)/r;
%             end
%                  Source_half(ss,ff,FrameNo)=(A(:,ss)'*P_half(:,ff,FrameNo))/SorNum;
%             end
%         end
%         for ss = 1:SorNum
%             P(ss,:)=[Source_half(ss,:,FrameNo),zeros(1,1),fliplr(conj(Source_half(ss,2:end,FrameNo)))];
%             p_part(ss,:)=(ifft(P(ss,:),NFFT));
% 
%             % --overlap and add--
%             tt2 = 1:NWIN;
%             p(ss,t_start+tt2)=p(ss,t_start+tt2)+p_part(ss,tt2);
%         end
%     end

    
end