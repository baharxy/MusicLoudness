function [ logf , avlogA , pinkness , alpha , err , fiteval ] = Voss( a)
%VOSS Generate PSD of audio f i l e
% [logf, logA, alpha, err, fit] = Voss( filename, minutes )
% Generates the power spectral density and value of alpha for an input
% audio file using the instentaneous loudness.


% Converts to a mono single and squares
Fs=1/(.01);
A= abs (fft (a)) ;
% Generates fourier transform
as=numel(a);
A=A ( 1 : round (as/2) );
A=A .* conj (A) / as;
f =(0:round(as/2)-1)*(Fs/as); % Calculates frequencies
nA=sum(f<Fs/2);
f=f (2:nA);
A=A ( 2 : nA );
logf=log10 ( f );
logA=(log10(A))';
logf1=linspace (min(logf) ,max(logf) ,1000000);
logA=interp1(logf ,logA,logf1);
logf=logf1;
% Resamples data linearly
[linfit ,fitS]=polyfit(logf ,logA,1);
% Fits a linear line to the data
err=sqrt ( diag ( inv ( fitS .R)*inv (fitS .R') ) .* fitS . normr .^2./ fitS . df ); % Calculates the error on the gradient
alpha=linfit(1);
% Saves the gradient as alpha
err=err(1);
fiteval=polyval(linfit ,logf);
% Makes an array representing the fit
avlogA=smooth(logA,1000, 'lowess'); % Smooths the data
f=10.^logf;
oof=log10 (1./f ); 
oof=oof+(avlogA(1)-oof(1) );
cor=corrcoef (avlogA , oof );
cor=abs ( cor (1,2) );
if alpha<= -1; x=exp(alpha)/exp(-1);
pinkness=cor *100*x;
else
x=exp(-alpha-2)/exp(-1); 
pinkness=cor *100*x;
end
% Calculates the pinkness of the data , a number taking the corrolation and
% the gradient into account
end

    