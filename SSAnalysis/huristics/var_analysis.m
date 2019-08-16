function [H ,x,y, yfit]=var_analysis(sequence)

% long_v=length(X);
% start=log10(2);
% stop=floor(log10((long_v/2)));
% number_of_points=floor(long_v/20);
% lag=floor(logspace(start,stop, number_of_points));
% result1=[]; size=[];
% for i=1:length(lag)
%     t_subgr=lag(i);
%     k=floor(long_v/t_subgr); %number of blocks
%     result=[];
%     for j=1:k
%         first=(j-1)*(t_subgr)+1;
%         x1=mean(X(first:first+t_subgr-1));
%         result=[result x1];
%     end
%     result1=[result1 var(result)];
%     size=[size t_subgr];
% end
% inds=find(diff(log10(result1'))>0.01);
% % figure(); plot(log10(size), log10(result1), 'b.'); hold on;
%  p= polyfit(log10(size(inds)), log10(result1(inds)),1);
% % f=polyval(p, log10(size(inds)));
% %plot(log10(size(inds)), f, 'r--')
% H=p(1)/2+1;


N = length(sequence);
mlarge = floor(N/5);
M = [floor(logspace(0,log10(mlarge),50))];
M = unique(M(M>1));
n = length(M);
cut_min = ceil(n/10);
cut_max = floor(6*n/10);
V = zeros(1,n);
for i = 1:n
    m = M(i);
    k = floor(N/m);
    matrix_sequence = reshape(sequence(1:m*k),m,k);
    V(i) = var(sum(matrix_sequence,1)/m);
end
x = log10(M);
y = log10(V);
y1 = -x+y(1)+x(1);
X = x(cut_min:cut_max);
Y = y(cut_min:cut_max);
p1 = polyfit(X,Y,1);
Yfit = polyval(p1,X);
yfit = polyval(p1,x);
alpha = -(Yfit(end)-Yfit(1))/(X(end)-X(1));
H = 1-alpha/2;