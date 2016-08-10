%2-phase Monte Carlo light escape program
%Script that runs a Monte Carlo approach to find the emission spectrum of a
%2-species mix, given the single species absorption and emission spectra. 

tic
%%_____Variables to play with
nphot=5e5;              %number of photons considered
dx=3e-3/2; dy=3e-2/2;   %size of the cuvette, x and y dimensions. Fate of photons leaving in the y direction
                            %may want to be changed. 
PLQEA=0.34;             %PLQE of substance A, currently PbCsCl3
PLQEB=0.97;             %PLQE of substance B, PbCsI3  
concA= 0.021*50;           %Concentrations of A {B}, in mg/ml
concB= 0.021;
conversionFactor=0.1;   %Converts the stored absorption constant into the same units as we are using in the model

excitationRatio = 0.117;   %P(initially excite A)/P(excite B) - the ratio of absorptions at the excitation wavelength

%The absorption length Abs = absorbance/length/log_10(e) and PL=PL.dL are stored
%in files, as is the wavelngth index wvI. 
load('PLandAbsA_PbCsCl3_B_PbCsI3.mat');
%Cut down range of PL/absorption data?
Lwwv=111;       %  lower cutoff of wavelength. 1 if all data wanted. 
Upwv = 700;     %  upper cutoff of data. length(wvI) if all wanted. 

%%____Initialisation

wvI=wvI(Lwwv:Upwv);
PLA=PLA(Lwwv:Upwv);
PLA=PLA/sum(PLA);
PLB=PLB(Lwwv:Upwv);
PLB=PLB/sum(PLB);
%Also convert the stored absorption into alphas/g into actual alphas
AlphaA=AbsA(Lwwv:Upwv)*concA*conversionFactor;
AlphaB=AbsB(Lwwv:Upwv)*concB*conversionFactor;
theta=acos(rand(nphot,1));  %random angles
esc=nan*ones(nphot,1);      %have we escaped (yet)?
xs=zeros(nphot, 1);         %Where are we?
ys=zeros(nphot, 1);
wvs=[ProbWv(round(excitationRatio*concA*nphot*PLQEA/(excitationRatio*concA*PLQEA+concB*PLQEB)),PLA);...
     ProbWv(round(concB*nphot*PLQEB/(excitationRatio*concA*PLQEA+concB*PLQEB)),PLB)];
AbsMed=false(nphot,1);      %Initialise the absorption medium recorder - 0 for A, 1 for B. 

%%____Calculation: 
%given photons travelling in a given direction, how far
%do they go and what substance absorbs them first? 
while any(isnan(esc))
    ins=isnan(esc);     %Selects photons still in play
    [xs(ins),ys(ins),AbsMed(ins)]=FateOfPhoton(...
        wvs(ins),theta(ins),AlphaA,AlphaB,xs(ins),ys(ins));    
    outs=(abs(xs)>dx);
    esc(outs)=1;            %Light that is detected
    esc(abs(ys)>dy)=0;      %light leaving the other side is lost
    ins=isnan(esc);
    Aswitch=AbsMed(ins);
    Bswitch=(AbsMed(ins)==0);
    wvs(Aswitch)=(rand(sum(Aswitch),1)<PLQEA).*ProbWv(sum(Aswitch),PLA);
    wvs(Bswitch)=(rand(sum(Bswitch),1)<PLQEB).*ProbWv(sum(Bswitch),PLB);
    theta=acos(2*rand(nphot,1)-1); 
    esc(wvs==0)=0;
end
StatResults=zeros(length(wvI),1);       %These are the final, wavelength-dependent results
for ii=1:length(wvI)
    StatResults(ii)=sum(wvs==ii);
end
plot(wvI,StatResults/max(StatResults))
xlabel('Wavelength (nm)')
ylabel('Photoluminescence (A.U.)')
toc