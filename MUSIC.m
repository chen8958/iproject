function [angel_curve]=MUSIC(MicPos,SorPos)
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

Rxx=zeros(MicNum,MicNum,length(Freqs));


for FrameNo=1:20
t_start=(FrameNo-1)*hopsize;
tt=(t_start+1):(t_start+NWIN);
for i=1:MicNum
    p_fft(i,tt)=fft(p_source(i,tt));
end
for ff=1:1:length(Freqs)
    %Rxx_tmp(:,:,ff) = p_fft(:,1:512)*p_fft(:,1:512)';
    Rxx_tmp(:,:,ff) = p_fft(:,:)*p_fft(:,:)';
end
Rxx=Rxx+Rxx_tmp;
end
Ps=zeros(MicNum,MicNum,length(Freqs));
P=zeros(MicNum,MicNum,length(Freqs));
for ff=1:1:length(Freqs);
[V,D]=eig(Rxx(:,:,ff));
[D_sort,index] = sort(diag(D),'descend');
% D_sort = D_sort(index);
V_sort = V(:,index);
%display(['now = ' num2str(ff)]);
%display(D_sort);
% for SignalSub=1:3
% P(:,:,ff)=P(:,:,ff)+V_sort(:,SignalSub)*V_sort(:,SignalSub)';
% end
%noise
%後面是noise subspace
for SignalSub=1:14
%v=V_sort(:,SignalSub)*l*V_sort(:,SignalSub)/(norm(l*V_sort(:,SignalSub))^2);
%P(:,:,ff)=P(:,:,ff)+v*v';

%P(:,:,ff)=P(:,:,ff)+(eye(MicNum)-V_sort(:,SignalSub)*V_sort(:,SignalSub)');
Ps(:,:,ff)=Ps(:,:,ff)+(V_sort(:,SignalSub)*V_sort(:,SignalSub)');
end
P(:,:,ff)=eye(MicNum)-Ps(:,:,ff);
end

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
    angel_curve=zeros(90,360);
    %angel_plant=zeros(length(Freqs),360);
    
    
    % 原music
%     for deg=1:360
%         for ff=1:length(Freqs)
%             
%             a_real=zeros(1,MicNum);
%             for ss=1:SorNum
%             k = 2*pi*Freqs(ff)/c; 
%             kappa = [cosd(SorPos(ss,1)).*sind(SorPos(ss,2)) sind(SorPos(ss,1)).*sind(SorPos(ss,2)) cosd(SorPos(ss,2))];
%             for MicNo=1:MicNum
%                 a(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
%             end
%                 a_real=a_real+a;
%             end
%                 
%             k = 2*pi*Freqs(ff)/c; 
%             kappa=[cosd(deg),sind(deg),0];
%             for MicNo=1:MicNum
%                 a(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
%             end
%              w=(inv(Rxx(:,:,ff)+0.01*eye(6))*a.')/(conj(a)*inv(P(:,:,ff)+0.01*eye(6))*a.');
% %             w=(inv(Rxx(:,:,ff)+0.01*eye(6))*a.')/(conj(a)*inv(Rxx(:,:,ff)+0.01*eye(6))*a.');
%             
%             angel_curve(deg)=angel_curve(deg)+abs(w'*a_real.');
%             angel_plant(ff,deg)=angel_plant(ff,deg)+abs(w'*a_real.');
%         end
%     end

    
    for ff=1:length(Freqs)   
        for ele=1:90
        for azi=1:360    
%             a_real=zeros(1,MicNum);
%             for ss=1:SorNum
%             k = 2*pi*Freqs(ff)/c; 
%             kappa = [cosd(SorPos(ss,1)).*sind(SorPos(ss,2)) sind(SorPos(ss,1)).*sind(SorPos(ss,2)) cosd(SorPos(ss,2))];
%             for MicNo=1:MicNum
%                 a(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
%             end
%                 a_real=a_real+a;
%             end
                
            k = 2*pi*Freqs(ff)/c; 
            kappa=[cosd(azi)*cosd(ele),sind(azi)*cosd(ele),sind(ele)];
            a=zeros(MicNum,1);
            for MicNo=1:MicNum
                a(MicNo)=exp(1i*k*kappa*MicPos(:,MicNo));
            end
             
             w=1/(a'*P(:,:,ff)*a);
%             w=1/(a'*inv(0.1*eye(6)-P(:,:,ff))*a);
%             w=(inv(Rxx(:,:,ff)+0.01*eye(6))*a.')/(conj(a)*inv(Rxx(:,:,ff)+0.01*eye(6))*a.');
            
            angel_curve(ele,azi)=angel_curve(ele,azi)+abs(w);
            %angel_plant(ff,deg)=angel_plant(ff,deg)+abs(w);
            %display([' azimuth = ' num2str(azi) ' elevation = ' num2str(ele) ' freq = ' num2str(ff) '/' num2str(length(Freqs))]);
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
    %figure(2)
    %plot(deg,angel_curve);
    
    

end