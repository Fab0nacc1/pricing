function [price,table] = barrier_fd(s0,X,T,hu,hd,r,q,sigma,M,n)
tic
dt = T/M;
bb = r-q;
dp = (hu-hd)/n;
table = zeros(M+1,n+1);
table(M+1,2:n ) = max(hd+dp*(1:n-1)-X,0);
a = dt./(2*(1+r.*dt)).*(sigma^2.*(0:n).^2-bb.*(0:n));
b = dt./(1+r.*dt).*(1./dt-sigma^2.*(0:n).^2);
c = dt./(2*(1+r.*dt)).*(sigma^2.*(0:n).^2+bb.*(0:n));
A = diag(b) + diag(a(2:n+1), -1) + diag(c(1:n), 1);
A = A(2:n,:);
% disp([dt,dp]);
% disp(A);
for j=M:-1:1
    %disp(table(j+1,:)');
    v = A*table(j+1,:)';
    table(j,2:n) = v';
end
sub = int32(floor((s0-hd)./dp));
w = (s0-hd)./dp-double(sub)*1.0;
price = (1-w)*table(1,sub+1)+w*table(1,sub+2);
toc
    

