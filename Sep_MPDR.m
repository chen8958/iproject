function Sep_MPDR(SorNum,MicPos,SorPos)
    fs=16000;
    c=343.0;
    %kappa=[sind(SorPos(:,1)).*sind(SorPos(:,2)),cosd(SorPos(:,1)).*sind(SorPos(:,2)),cosd(SorPos(:,2))];
    %kappa = [cosd(SorPos(1))*sind(SorPos(2)) sind(SorPos(1))*sind(SorPos(2)) cosd(SorPos(2))];
    kappa = [cosd(SorPos(:,1)).*cosd(SorPos(:,2)) sind(SorPos(:,1)).*cosd(SorPos(:,2)) sind(SorPos(:,2))];
    [D MicNum]=size(MicPos);
    for i=1:MicNum
        [p_source(i,:) fs]=audioread("p"+i+".wav");
    end
    
     %% Windowing
    NWIN=1024;
    hopsize=NWIN/2;                                                            % 50% overlap
    NumOfFrame=2*floor(length(p_source)/NWIN)-1;                                         % number of frames
    win = hann(NWIN+1);                                                        % hanning window
    win = win(1:end-1).';

    %% FFT
    NFFT=2^nextpow2(NWIN);
    df=fs/NFFT;
    Freqs=0:df:(NFFT/2-1)*df;
%   %% find Rxx
%     for i=1:MicNum
%         p_fft(i,:)=fft(p_source(i,:));
%     end
%     Rxx = p_fft*p_fft'+0.1*eye(6);
%%
%MVDR
% mu=1*pi;           %for request 1 pi 2*pi
% for m=1:6
%     for n=1:6
%         Rxx(m,n)=sinc((m-n)*mu/pi);
%     end
% end

Rxx=zeros(MicNum,MicNum,length(Freqs));

for FrameNo=1:20
t_start=(FrameNo-1)*hopsize;
tt=(t_start+1):(t_start+NWIN);
for i=1:MicNum
    p_fft(i,tt)=fft(p_source(i,tt));
end
for ff=1:1:length(Freqs)
    Rxx_tmp(:,:,ff) = p_fft(:,1:512)*p_fft(:,1:512)';
end
Rxx=Rxx+Rxx_tmp;
end

    %%    
    p=zeros(SorNum,length(p_source));
    for MicNo=1:MicNum
        for FrameNo=1:NumOfFrame
            t_start=(FrameNo-1)*hopsize;
            tt=(t_start+1):(t_start+NWIN);
            %P_half(MicNo,:,FrameNo)=fft(p_source(MicNo,tt).*win);
            P_whole(MicNo,:,FrameNo)=fft(p_source(MicNo,tt).*win);
            P_half(MicNo,:,FrameNo)=P_whole(MicNo,1:1:NFFT/2,FrameNo);
        end
    end
    for FrameNo=1:NumOfFrame
        % --time segment--
        t_start=(FrameNo-1)*hopsize;
        tt=(t_start+1):(t_start+NWIN);
          %% find Rxx
%       for i=1:MicNum
%           p_fft(i,:)=fft(p(i,tt));
%       end
%       Rxx = p_fft*p_fft'+0.01*eye(6);
        % --propagation--
        for ff=1:length(Freqs)
            k = 2*pi*Freqs(ff)/c;                                              
            for ss = 1:SorNum
            for m = 1:MicNum
                 
                 A(m,ss) =exp(1j*k*kappa(ss,:)*MicPos(:,m));

            end
                 %w=(inv(Rxx(:,:,ff)+0.001*eye(size(A,1)))*A(:,ss))/(A(:,ss)'*inv(Rxx(:,:,ff)+0.001*eye(size(A,1)))*A(:,ss));
                 w=inv(Rxx(:,:,ff)+0.01*eye(MicNum))*A(:,ss)/(A(:,ss)'*inv(Rxx(:,:,ff)+0.01*eye(MicNum))*A(:,ss));
                 Source_half(ss,ff,FrameNo)=(w'*P_half(:,ff,FrameNo));
%                     Source_half(ss,ff,FrameNo)=(A(:,ss)'*P_half(:,ff,FrameNo))/SorNum;
            end
        end
        for ss = 1:SorNum
            P(ss,:)=[Source_half(ss,:,FrameNo),zeros(1,1),fliplr(conj(Source_half(ss,2:end,FrameNo)))];
%             P(ss,:)=[Source_half(ss,:,FrameNo)];
            p_part(ss,:)=(ifft(P(ss,:),NFFT));

            % --overlap and add--
            tt2 = 1:NWIN;
            p(ss,t_start+tt2)=p(ss,t_start+tt2)+p_part(ss,tt2);
        end
    end
    for ss = 1:SorNum
        audiowrite(['MPDR_sep' num2str(MicNum) num2str(ss) '.wav'],real(p(ss,:))/max(abs(p(ss,:))),fs);
    end
end