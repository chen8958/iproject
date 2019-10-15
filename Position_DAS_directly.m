function [angel_plant]=Position_DAS_directly(MicPos,SorPos)
c=343;
path=cd;
SorNum=size(SorPos,1);
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

%     for i=1:MicNum
%     x(i,1:4*fs)=p(i,1:4*fs);
%     fft_x(i,:)=fft(x(i,:),NFFT);
%     end
    
    
    angel_curve=zeros(1,360);
    angel_plant=zeros(length(Freqs),360);
    for deg=1:360
        for ff=1:length(Freqs)
%             k = 2*pi*Freqs(ff)/c; 
%             kappa = [cosd(SorPos(1))*sind(SorPos(2)) sind(SorPos(1))*sind(SorPos(2)) cosd(SorPos(2))];
%             for MicNo=1:MicNum
%                 a_real(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
%             end
            

%         k = 2*pi*Freqs(ff)/c; 
%         kappa1 = [cosd(SorPos(1,1)).*sind(SorPos(1,2)) sind(SorPos(1,1)).*sind(SorPos(1,2)) cosd(SorPos(1,2))];
% %         kappa1 = [cosd(-30) sind(-30) 0];     % look direction
%         kappa2 = [cosd(SorPos(2,1)).*sind(SorPos(2,2)) sind(SorPos(2,1)).*sind(SorPos(2,2)) cosd(SorPos(2,2))];
% %         kappa2 = [cosd(150) sind(150) 0];
%         a1 =exp(1j*k*kappa1*MicPos).';                           % manifold vector
%         a2 =exp(1j*k*kappa2*MicPos).';
%         a_real=a1+a2;


%             k = 2*pi*Freqs(ff)/c; 
%             kappa1 = [cosd(SorPos(1,1)).*sind(SorPos(1,2)) sind(SorPos(1,1)).*sind(SorPos(1,2)) cosd(SorPos(1,2))];
%             a1=exp(1i*k*kappa1*MicPos);
%             kappa2 = [cosd(SorPos(2,1)).*sind(SorPos(2,2)) sind(SorPos(2,1)).*sind(SorPos(2,2)) cosd(SorPos(2,2))];
%             a2=exp(1i*k*kappa2*MicPos);
%             a_real=a1+a2;
            a_real=zeros(1,MicNum);
            for ss=1:SorNum
            k = 2*pi*Freqs(ff)/c; 
            kappa = [cosd(SorPos(ss,1)).*sind(SorPos(ss,2)) sind(SorPos(ss,1)).*sind(SorPos(ss,2)) cosd(SorPos(ss,2))];
            for MicNo=1:MicNum
                a_tmp(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
            end
                a_real=a_real+a_tmp;
            end
            
            k = 2*pi*Freqs(ff)/c; 
            kappa=[cosd(deg),sind(deg),0];
            for MicNo=1:MicNum
                a(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
            end
            w=conj(a);
            angel_curve(deg)=angel_curve(deg)+abs(w*a_real.');
            angel_plant(ff,deg)=abs(w*a_real.');

% 
%         k = 2*pi*Freqs(ff)/c; 
% %         kappa1 = [sind(-30) cosd(-30) 0];     % look direction
% %         kappa2 = [sind(150) cosd(150) 0];
%         
%         kappa1 = [cosd(-30) sind(-30) 0];     % look direction
%         %kappa2 = [cosd(150) sind(150) 0];
%         a =exp(1j*k*kappa1*MicPos).';                           % manifold vector
%         %a2 =exp(1j*k*kappa2*MicPos).';
%         %a=a1+a2;
%         theta =[cosd(deg),sind(deg),0];                              % source direction  
%         W= exp(1i*k*theta*MicPos).';
%         angel_plant(ff,deg) = abs(W'*a);
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