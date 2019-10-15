function [P_half SorPos SorLen p]=Mix3D_Plan_function_noise(MicPos,SorPos)
c=343;  % sound speed
[N M]=size(MicPos);
[SorNum D]=size(SorPos);
path=cd;
% MicPos=MicPos+(rand(N,M)-rand(N,M))*0.3;

%% Sound source(s)
cd([path '\audio_R'])
[x1 fs]=audioread('female_16k_10s.wav');
[x2 fs]=audioread('male_16k_10s.wav');
% y1=wgn(length(x1),1,10);
% x1=x1+y1;
% y2=wgn(length(x2),1,10);
% x2=x2+y2;

% x1=awgn(x1,30);
% x2=awgn(x2,30);

% [x3 fs]=audioread('4_speech_5s.wav');
cd(path)
SorLen=fs*4;
%Source=[x1(1:fs*4)];
Source=[x1(1:fs*4),x2(1:fs*4)];
% Source=[x1(1:fs*1),x2(1:fs*1),x3(1:fs*1)];

%% Add noise to each microphone
for MicNo=1:M
    for ss=1:SorNum
%         source(ss,:,MicNo)=awgn(Source(:,ss)',40);
        source(ss,:,MicNo)=Source(:,ss)';
    end
end

%% Windowing
NWIN=1024;
hopsize=NWIN/2;                                                            % 50% overlap
NumOfFrame=2*floor(4*fs/NWIN)-1;                                           % number of frames
win = hann(NWIN+1);                                                        % hanning window
win = win(1:end-1).';

%% FFT
NFFT=2^nextpow2(NWIN);
df=fs/NFFT;
Freqs=0:df:(NFFT/2-1)*df;

%%
p=zeros(M,fs*4);
for FrameNo=1:NumOfFrame
    % --time segment--
    t_start=(FrameNo-1)*hopsize;
    tt=(t_start+1):(t_start+NWIN);
    
    % --transform to frequency domain--
    for MicNo=1:M
        for ss=1:SorNum
            source_win(ss,:)=source(ss,tt,MicNo).*win;
            source_zp(ss,:)=[source_win(ss,:) zeros(1,(NFFT-NWIN))];
            SOURCE(ss,:)=fft(source_zp(ss,:),NFFT);
            SOURCE_half(ss,:,MicNo)=SOURCE(ss,1:NFFT/2);
        end
    end

    % --propagation--
    for ff=1:length(Freqs)
        k = 2*pi*Freqs(ff)/c;                                              
        for ss = 1:SorNum
            kappa = [cosd(SorPos(ss,1))*sind(SorPos(ss,2)) sind(SorPos(ss,1))*sind(SorPos(ss,2)) cosd(SorPos(ss,2))];
            A(:,ss)=exp(1i*k*kappa*MicPos).';
        end
        for MicNo=1:M
            P_half(MicNo,ff,FrameNo)=(A(MicNo,:)*SOURCE_half(:,ff,MicNo));   
        end
    end   
            
    for MicNo=1:M     
        P(MicNo,:)=[P_half(MicNo,:,FrameNo),zeros(1,1),fliplr(conj(P_half(MicNo,2:end,FrameNo)))];
        p_part(MicNo,:)=(ifft(P(MicNo,:),NFFT));
        
        %--overlap and add--
        p(MicNo,tt)=p(MicNo,tt)+p_part(MicNo,:);
    end
end

%% Write wave file
    %audiowrite(['p' num2str(1) '.wav'],p(1,:),fs);
for MicNo=1:M
    p(MicNo,:)=awgn(p(MicNo,:),30);
    audiowrite(['p' num2str(MicNo) '.wav'],p(MicNo,:)/max(abs(p(MicNo,:))),fs);
end
end

