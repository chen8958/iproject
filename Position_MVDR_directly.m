function [angel_plant]=Position_MVDR_directly(MicPos,SorPos)
c=343;
path=cd;
[D MicNum]=size(MicPos);
    %% Windowing
cd([path '\audio_R'])
[x1 fs]=audioread('female_16k_10s.wav');
SorLen=fs*4;
Source=[x1(1:fs*4)];
cd(path)
NWIN=1024;
hopsize=NWIN/2;                                                            % 50% overlap
NumOfFrame=2*floor(fs*4/NWIN)-1;                                           % number of frames
win = hann(NWIN+1);                                                        % hanning window
win = win(1:end-1).';
%% FFT
NFFT=2^nextpow2(NWIN);
df=fs/NFFT;
Freqs=0:df:(NFFT/2-1)*df;

%     for i=1:MicNum
%     x(i,1:4*fs)=p(i,1:4*fs);
%     fft_x(i,:)=fft(x(i,:),NFFT);
%     end
    
    
%     angel_curve=zeros(1,180);
%     angel_plant=zeros(length(Freqs),180);
%     for deg=1:180
%         for ff=1:length(Freqs)
%             k = 2*pi*Freqs(ff)/c; 
%             kappa = [cosd(SorPos(1))*cosd(SorPos(2)) sind(SorPos(1))*cosd(SorPos(2)) sind(SorPos(2))];
%             for MicNo=1:MicNum
%                 a_real(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
%             end
%             
%             k = 2*pi*Freqs(ff)/c; 
%             kappa=[cosd(deg),sind(deg),0];
%             for MicNo=1:MicNum
%                 a(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
%             end
%             w=conj(a);
%             %angel_curve(deg+90)=angel_curve(deg+90)+abs(w*fft_x(:,ff));
%             %angel_plant(ff,deg+90)=abs(w*fft_x(:,ff));
%             angel_curve(deg)=angel_curve(deg)+abs(w*a_real.');
%             angel_plant(ff,deg)=abs(w*a_real.');
%         end
%     end
%     figure(1)
%     deg=1:180;
%     pcolor(deg,Freqs,angel_plant);
%     shading interp;
%     figure(2)
%     plot(deg,angel_curve);
    
    
    
    
    angel_curve=zeros(90,360);
    for ff=1:length(Freqs)   
        for ele=1:90
        for azi=1:360    
            k = 2*pi*Freqs(ff)/c; 
            kappa = [cosd(SorPos(1))*cosd(SorPos(2)) sind(SorPos(1))*cosd(SorPos(2)) sind(SorPos(2))];
            for MicNo=1:MicNum
                a_real(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
            end
                
            k = 2*pi*Freqs(ff)/c; 
            kappa=[cosd(azi)*cosd(ele),sind(azi)*cosd(ele),sind(ele)];
            a=zeros(MicNum,1);
            for MicNo=1:MicNum
                a(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
            end
             
            w=1/(a'*P(:,:,ff)*a);
            angel_curve(ele,azi)=angel_curve(ele,azi)+abs(w);
           
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
    title('MUSIC');

end