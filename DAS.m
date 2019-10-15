%% shotgun model design by Delay And Sum (DAS)
close all;clear all;clc
%% parameters setting
fs=44100;                                                      % sampling rate
c=343;                                                         % sound speed
M=8;                                                           % number of microphone
d=0.01;                                                        % microphone interval
NFFT=1024;                                                     % cut 1024 frames
df=fs/NFFT;
Freqs=0:df:fs/2;                                               % effective frequency
EstiAng=-180:1:180;                                              % scan angle
r=0.015;
%% microphone position (ULA)
inter=360/M;
for i=1:M
   MicPos(:,i)=[r*cosd(inter*i); r*sind(inter*i)];
end
%% delay and sum algorithm
for ff=1:length(Freqs)                                         % scan frequency
    w=2*pi*Freqs(ff);
    k=w/c;
    for ang=1:length(EstiAng)                                  % scan angle
%         kappa = [sind(EstiAng(ang)) , cosd(EstiAng(ang))];     % look direction
%         a =exp(1j*k*kappa*MicPos).';                           % manifold vector
%         theta =[sind(0) cosd(0)];                              % source direction  
%         W= exp(1i*k*theta*MicPos).';
        kappa1 = [sind(-60) cosd(-60)];     % look direction
        kappa2 = [sind(120) cosd(120)];
        a1 =exp(1j*k*kappa1*MicPos).';                           % manifold vector
        a2 =exp(1j*k*kappa2*MicPos).';
        a=a1+a2;
        theta =[sind(EstiAng(ang)) , cosd(EstiAng(ang))];                              % source direction  
        W= exp(1i*k*theta*MicPos).';
        y(ff,ang) = abs(W'*a);
    end  
end
%% shotgun information
L=MicPos(1,8);
minlength=MicPos(1,2);
MicPos
%% diagram
for ff=1:length(Freqs)
    y_n(ff,:) = y(ff,:)./max(y(ff,:));    
end
figure(1)
plot(EstiAng, mean(y(ff,:),1)./max( mean(y(ff,:),1) ));
figure(2)
contourf(EstiAng,Freqs,y_n)
xlabel('Angle(deg)')
ylabel('Frequency(Hz)')
















