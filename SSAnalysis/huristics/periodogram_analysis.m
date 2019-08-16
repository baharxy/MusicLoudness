function [F, Pxx]=periodogram_analysis(X)
fs=1/(.01);
long_v=length(X);
[Pxx,  F]=periodogram(X,[],[], fs);
% figure(); plot(log10(F),log10(Pxx),'b.');
% polyfit(log10(F),log10(Pxx),1)