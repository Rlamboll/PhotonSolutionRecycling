function [ xout,yout,AbsMed] = FateOfPhoton( wv,theta,AbsA,AbsB,xin,yin)
%FATEOFPHOTON calculates how far a photon of wavelength wv and travelling at angle 
%theta travels in a medium with two
%mixed phases with given PL and absorptions. If |xin+dx | >d the photon
%escapes. The output is the new x, if it escaped and the absorbing medium. 

%INPUT : wavelengths, angles, absorption of the two species and starting x/y
%OUTPUT: end x,y and which medium absorbed. 

%Calculate length going before absorbed by particular materials
LA=-1./AbsA(wv).*log(rand(size(wv)));
LB=-1./AbsB(wv).*log(rand(size(wv)));

%Shorter length is the absorption material - A for true
AbsMed=LA<LB; 

xout(AbsMed)=xin(AbsMed)+cos(theta(AbsMed)).*LA(AbsMed);
xout(~AbsMed)=xin(~AbsMed)+cos(theta(~AbsMed)).*LB(~AbsMed);
yout(AbsMed)=yin(AbsMed)+sin(theta(AbsMed)).*(LA(AbsMed));
yout(~AbsMed)=yin(~AbsMed)+sin(theta(~AbsMed)).*(LB(~AbsMed));
end

