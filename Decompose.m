function [A,u2]=Decompose(u,v1,v2,wv0,wv1,wv2)
%Attempts to find the best fit to express one vector, u = A(1)*v1+A(2)*v2, where 
%u's values are defined over wv0, v1 has values at wv1 and v2 has values at
%wv2. 

%Returns A, the components of vectors v1 and v2 found in vector u for non-
%orthogonal v1,v2, and u2, the reconstructed best-fit sum. 

%Program should not need altering. 
    fitfn=fittype(@(B1,B2,x)B1.*interp1(wv1,v1,x,'linear',0)+B2*interp1(wv2,v2,x,'linear',0));
    fitResult=fit(wv0,u,fitfn,'StartPoint',[1,1]);
    u2=fitfn(fitResult.B1,fitResult.B2,wv0);
    A= [fitResult.B1 fitResult.B2];
end