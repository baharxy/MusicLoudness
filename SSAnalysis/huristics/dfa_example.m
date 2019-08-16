Y = randn(5000,1);
%X=Y;
X = %pinknoise(5000);% cumsum(Y);
%X=res_a.STL; %pinknoise(59232);%res_a.InstantaneousLoudness;
%Xp=res_c.STL;
%Xb=res_b.STL;
plot_fun = @(xp,A,ord) polyval(A,log(xp));
%% Example 1. DFA of order 1.
% A(1) is approximately 1.5. Indeed X is a brownian motion.
pts = 10:10:100;
[Aa,Fa] = DFA_fun(X,pts,1);
%[Ac,Fc] = DFA_fun(Xp,pts,2);
%[Ab,Fb] = DFA_fun(Xb,pts,2);
Aa(1)
figure()
scatter(log(pts),log(Fa)); hold on;
%scatter(log(pts),log(Fa)); hold on;
%scatter(log(pts),log(Fc))

hold on
x = 1:10:1000;
plot(log(x),plot_fun(x,Aa),'--')
hold off
%% Example 2. DFA of order 2
% A(1) is approximately 1.5. Indeed X is a brownian motion.
pts = 50:20:500;
[A,F] = DFA_fun(X,pts,2);
A(1)
figure()
scatter(log(pts),log(F))
hold on
x = 1:10:1000;
plot(log(x),plot_fun(x,A),'--')
hold off
%% Example 2. DFA of order 1
% A(1) is approximately 0.5. Indeed Y is an uncorrelated process.
pts = 10:10:200;
[A,F] = DFA_fun(Y,pts);
A(1)
figure()
scatter(log(pts),log(F))
hold on
x = 1:10:1000;
plot(log(x),plot_fun(x,A),'--')
hold off
%%