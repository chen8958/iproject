clc
clear all
%%
%for origin uca
%MicPos=(1/100)*[0,3/sqrt(2),3/sqrt(2),0,-3/sqrt(2),-3/sqrt(2);3,3/sqrt(2),-3/sqrt(2),-3,-3/sqrt(2),3/sqrt(2);0,0,0,0,0,0];
%%
%for origin uca
%MicPos=(1/100)*[4.5*cosd(120),4.5*cosd(60),4.5,4.5*cosd(-60),4.5*cosd(-120),-4.5;4.5*sind(120),4.5*sind(60),0,4.5*sind(-60),4.5*sind(-120),0;0,0,0,0,0,0];

%%
% % for 24 mic
% angel=0:45:315;
% MicPos=(1/100)*4.5*[cosd(angel);sind(angel);zeros(1,length(angel))]

%%

% angel=0:45:315;
% MicPos=(1/100)*2.6*[cosd(angel);sind(angel);zeros(1,length(angel))]

%%
%uca+ula
% angel=0:60:300;
% MicPos=(1/100)*4.5*[cosd(angel);sind(angel);zeros(1,length(angel))]
% MicPos2=(1/100)*[0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;1,2,3,4,5,6,7,8];
% MicPos=[MicPos,MicPos2];

%%
%usa
% angel=0:60:300;
% MicPos=(1/100)*4.5*[cosd(angel);sind(angel);zeros(1,length(angel))]
% 
% angel2=0:72:288;
% elev2=60.*ones(1,length(angel2))
% MicPos2=(1/100)*4.5*[cosd(angel2).*sind(elev2);sind(angel2).*sind(elev2);cosd(elev2)];
% 
% angel3=0:180:180;
% elev3=30.*ones(1,length(angel3))
% MicPos3=(1/100)*4.5*[cosd(angel3).*sind(elev3);sind(angel3).*sind(elev3);cosd(elev3)];
% 
% MicPos4=(1/100)*4.5*[0;0;1];
% 
% MicPos=[MicPos,MicPos2,MicPos3,MicPos4];

%%
%ula 3cm
angel=0:180:180;
MicPos=(1/100)*1.5*[cosd(angel);sind(angel);zeros(1,length(angel))];


SorPos=[180,90];