function [angel_curve]=MVDR(MicPos,SorPos)
c=343;
path=cd;
SorNum=size(SorPos,1);
[D MicNum]=size(MicPos);
    %% Windowing
      for i=1:MicNum
          [p_source(i,:) fs]=audioread("p"+i+".wav");
      end
    %% Windowing
% cd([path '\audio_R'])
% [x1 fs]=audioread('female_16k_10s.wav');
% SorLen=fs*4;
% Source=[x1(1:fs*4)];
% cd(path)

%%
%windowing
NWIN=1024;

hopsize=NWIN/2;                                                            % 50% overlap
NumOfFrame=2*floor(fs*4/NWIN)-1;                                           % number of frames
win = hann(NWIN+1);                                                        % hanning window
win = win(1:end-1).';
%% FFT
NFFT=2^nextpow2(NWIN);
df=fs/NFFT;
Freqs=0:df:(NFFT/2-1)*df;
%%find Rxx



% 
% for i=1:MicNum
%     p_fft(i,:)=fft(p_source(i,:));
% end
% for ff=1:1:length(Freqs)
%     Rxx(:,:,ff) = p_fft(:,ff)*p_fft(:,ff)';
% end


%% MVDR
           %for request 1 pi 2*pi
% for m=1:6
%     for n=1:6
%         Rxx(m,n)=sinc((norm(MicPos(:,m)-MicPos(:,n)))/pi);
%     end
% end
% figure(3)
% mesh(Rxx)


%     for i=1:MicNum
%     x(i,1:4*fs)=p(i,1:4*fs);
%     fft_x(i,:)=fft(x(i,:),NFFT);
%     end
%     angel_curve=zeros(1,360);
%     angel_plant=zeros(length(Freqs),360);
%     for deg=1:360
%         for ff=1:length(Freqs)
% %             k = 2*pi*Freqs(ff)/c; 
% %             kappa = [cosd(SorPos(1,1)).*sind(SorPos(1,2)) sind(SorPos(1,1)).*sind(SorPos(1,2)) cosd(SorPos(1,2))];
% %             for MicNo=1:MicNum
% %                 a_real(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
% %             end
%             a_real=zeros(1,MicNum);
%             for ss=1:SorNum
%             k = 2*pi*Freqs(ff)/c; 
%             kappa = [cosd(SorPos(ss,1)).*sind(SorPos(ss,2)) sind(SorPos(ss,1)).*sind(SorPos(ss,2)) cosd(SorPos(ss,2))];
%             for MicNo=1:MicNum
%                 a_tmp(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
%             end
%                 a_real=a_real+a_tmp;
%             end
%             
%             k = 2*pi*Freqs(ff)/c; 
%             kappa=[cosd(deg),sind(deg),0];
%             for MicNo=1:MicNum
%                 a(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
%             end
%             
%             for m=1:MicNum
%                 for n=1:MicNum
%                     Rnn(m,n)=sinc((norm(MicPos(:,m)-MicPos(:,n)))*k/pi);
%                 end
%             end
%             
%             
%             %w=(inv(Rnn+0.01*eye(MicNum))*a.')/(conj(a)*inv(Rnn+0.01*eye(MicNum))*a.');
%             w=(inv(Rnn+0.01*eye(MicNum))*a.')/(conj(a)*inv(Rnn+0.01*eye(MicNum))*a.');
%             
%             %angel_curve(deg+90)=angel_curve(deg+90)+abs(w*fft_x(:,ff));
%             %angel_plant(ff,deg+90)=abs(w*fft_x(:,ff));
%             angel_curve(deg)=angel_curve(deg)+abs(w'*a_real.');
%             angel_plant(ff,deg)=abs(w'*a_real.');
%         end
%     end
%     figure(1)
%     deg=1:360;
%     pcolor(deg,Freqs,angel_plant);
%     shading interp;
%     figure(2)
%     plot(deg,angel_curve);

% for m=1:MicNum
% for n=1:MicNum
%     Rnn(m,n)=sinc((norm(MicPos(:,m)-MicPos(:,n)))*k/pi);
% end
% end


    
    angel_curve=zeros(90,360);
    for ff=1:length(Freqs)
        
        for s=1:SorNum
        k = 2*pi*Freqs(ff)/c; 
        kappa = [cosd(SorPos(s,1))*cosd(SorPos(s,2)) sind(SorPos(s,1))*cosd(SorPos(s,2)) sind(SorPos(s,2))];
        for MicNo=1:MicNum
            a_real(s,MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
        end
        end
        
        for m=1:MicNum
        for n=1:MicNum
             Rnn(m,n)=sinc((norm(MicPos(:,m)-MicPos(:,n)))*k/pi);
        end
        end

        
        for ele=1:90
        for azi=1:360    

            k = 2*pi*Freqs(ff)/c; 
            kappa=[cosd(azi)*cosd(ele),sind(azi)*cosd(ele),sind(ele)];
            a=zeros(MicNum,1);
            for MicNo=1:MicNum
                a(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
            end
            
            
             
%             w=1/(a'*P(:,:,ff)*a);
            w=(inv(Rnn+0.01*eye(MicNum))*a)/(a'*inv(Rnn+0.01*eye(MicNum))*a);
            
            for s=1:SorNum
            angel_curve(ele,azi)=angel_curve(ele,azi)+abs(w'*a_real(s,:).');
            end
%             angel_curve(ele,azi)=angel_curve(ele,azi)+abs(w'*a_real.');
           
        end
        end
        display([' freq = ' num2str(ff) '/' num2str(length(Freqs))]);
    end
    
    figure(1);
    
%     pcolor(deg,Freqs,angel_plant);
    %shading interp;
    azi=1:360;
    ele=1:90;
    contourf(azi,ele,angel_curve);
    
    hold on;
    for i=1:SorNum
    plot(SorPos(i,1),SorPos(i,2),'x','LineWidth',5,'MarkerSize',15,'Color',[1,0,0]);
    end
    
    xlabel('azi');
    ylabel('ele');
    title('MVDR');
    

end