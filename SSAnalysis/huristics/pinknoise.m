%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Pink Noise Generation with MATLAB Implementation   %
%                                                      %
% Author: M.Sc. Eng. Hristo Zhivomirov       07/30/13  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = pinknoise(N)
B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
A = [1 -2.494956002   2.017265875  -0.522189400];
nT60 = round(log(1000)/(1-max(abs(roots(A))))); % T60 est.
v = randn(1,N+nT60); % Gaussian white noise: N(0,1)
x = filter(B,A,v);    % Apply 1/F roll-off to PSD
y = x(nT60+1:end);    % Skip transient response
% function: y = pinknoise(N) 
% N - number of samples to be returned in row vector
% y - row vector of pink (flicker) noise samples

% The function generates a sequence of pink (flicker) noise samples. 
% Pink noise has equal energy in all octaves (or similar log bundles) of frequency.
% In terms of power at a constant bandwidth, pink noise falls off at 3 dB per octave. 
% http://www.mathworks.com/matlabcentral/fileexchange/42919-pink--red--blue-and-violet-noise-generation-with-matlab-implementation/content/example.m

% difine the length of the vector
% ensure that the M is even
% if rem(N,2)
%     M = N+1;
% else
%     M = N;
% end
% 
% % generate white noise
% x = randn(1, M);
% 
% % FFT
% X = fft(x);
% 
% % prepare a vector for 1/f multiplication
% NumUniquePts = M/2 + 1;
% n = 1:NumUniquePts;
% n = sqrt(n);
% 
% % multiplicate the left half of the spectrum so the power spectral density
% % is proportional to the frequency by factor 1/f, i.e. the
% % amplitudes are proportional to 1/sqrt(f)
% X(1:NumUniquePts) = X(1:NumUniquePts)./n;
% 
% % prepare a right half of the spectrum - a copy of the left one,
% % except the DC component and Nyquist frequency - they are unique
% X(NumUniquePts+1:M) = real(X(M/2:-1:2)) -1i*imag(X(M/2:-1:2));
% 
% % IFFT
% y = ifft(X);
% 
% % prepare output vector y
% y = real(y(1, 1:N));
% 
% % ensure unity standard deviation and zero mean value
% y = y - mean(y);
% yrms = sqrt(mean(y.^2));
% y = y/yrms;

end