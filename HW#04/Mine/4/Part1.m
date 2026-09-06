clc
clear
close all
c = [0.2 0.4 0.6 -0.1 -0.3];
d = [0.1 0.3 -0.2 0.5 -0.3];
fs = 20;
t = @(k)(k-1 : 1/fs : k -1/fs);
s1k = @(k)c(k)*sin(2*pi*t(k));
s2k = @(k)d(k)*sin(4*pi*t(k));
A = [0.8 -0.6;0.6 0.8];
s1 = [s1k(1),s1k(2),s1k(3),s1k(4),s1k(5)];
s2  = [s2k(1),s2k(2),s2k(3),s2k(4),s2k(5)];
S = [s1;s2];
X = A*S;
s1 = s1/norm(s1);
s2 = s2/norm(s2);
figure(1)
plot((0:length(s1)-1)/fs,s1)
title("source1 Signal time domain")
figure(2)
plot((0:length(s2)-1)/fs,s2)
title("source2 Signal time    domain")
figure(3)
plot((0:length(X(1,:))-1)/fs,X(1,:))
title("X1t Signal time domain")
figure(4)
plot((0:length(X(2,:))-1)/fs,X(2,:))
title("X2t Signal time domain")
hold off

                        %%%%
Rx = X*X';
[U,D] = eig(Rx);
Z = sqrt(inv(D))*U.'*X;
Z1 = Z(:,1:20);
Z2 = Z(:,21:40);
RZ1 = Z1*Z1.';
RZ2 = Z2*Z2.';
[Q,~] = eig(RZ1,RZ2);
B = Q.';
sHat = B*Z;

sHat(1,:) = sHat(1,:)/norm(sHat(1,:));
sHat(2,:) = sHat(2,:)/norm(sHat(2,:));

if norm(sHat(1,:)-s1) >= 0.5 && norm(-sHat(1,:)-s1) >= 0.5
sHat = [sHat(2,:);sHat(1,:)];
end
if norm(sHat(1,:)-s1) >= norm(-sHat(1,:)-s1)
    sHat(1,:) = -sHat(1,:);
end
if norm(sHat(2,:)-s2) >= norm(-sHat(2,:)-s2)
    sHat(2,:) = -sHat(2,:);
end
Error = sHat-[s1;s2];
normError = norm(Error);
ErrorPb = (normError / norm([s1;s2]))^2

 
              %%%%%
 ErrorPc = 0;
 K = 5;
 B = [1 0 ;0 0] ;
for itr = 1:100
 for i = 1 : 2
   Ri = zeros(2);
  for k = 1:K
     ind = setdiff(1:2   ,i);
     Zi = Z(:,20*k - 19:20*k);
     Rzi = Zi*Zi.';
     Ri = Ri +(Rzi*B(ind,:).'* B(ind,:) *Rzi.');
  end  
     [V,lambdai] = eig(Ri,'vector');
     [lambdai,indi] = sort(lambdai);
     V = V(:,indi);
     B(i,:) = V(:,1).';
   if i > 1
      Rb = B(1:i-1,:).'*B(1:i-1,:);
      B(i,:) = (eye(2)-Rb)*(B(i,:).');
      B(i,:) = B(i,:)/norm(B(i,:));
   end
 end
 sHat = B*Z;
sHat(1,:) = sHat(1,:)/norm(sHat(1,:));
sHat(2,:) = sHat(2,:)/norm(sHat(1,:));


if norm(sHat(1,:)-s1) >= 0.5 && norm(-sHat(1,:)-s1) >= 0.5
sHat = [sHat(2,:);sHat(1,:)];
end
if norm(sHat(1,:)-s1) >= norm(-sHat(1,:)-s1)
    sHat(1,:) = -sHat(1,:);
end
if norm(sHat(2,:)-s2) >= norm(-sHat(2,:)-s2)
    sHat(2,:) = -sHat(2,:);
end
Error = sHat-[s1;s2];
normError = norm(Error);
ErrorPc = ErrorPc + (normError/norm([s1;s2]))^2;
end
ErrorPc = ErrorPc/100

                     %%%
W = randn(size(X));
W = W/ norm(W);
sig = sqrt(norm(X)^2/100);
Y = X+sig*W;
[U,D] = eig(Y*Y.');
Z = sqrt(inv(D))*U.'*Y;
K =5;
B = [0 0 ; 0 0];
ErrorPd = 0;
for itr = 1:100
 for i = 1 : 2
   Ri = zeros(2);
  for k = 1:K
     ind = setdiff(1:2,i);
     Zi = Z(:,20*k - 19:20*k);
     Rzi = Zi*Zi.';
     Ri = Ri +(Rzi*B(ind,:).'* B(ind,:) *Rzi.');
  end  
     [V,lambdai] = eig(Ri);
     [lambdai,indi] = sort(lambdai);
     V = V(:,indi);   
     B(i,:) = V(:,1).';
     
  if i > 1
     Rb = B(1:i-1,:).'*B(1:i-1,:);
      B(i,:) = (eye(2)-Rb)*(B(i,:).');
      B(i,:) = B(i,:)/norm(B(i,:));
  end
 end
sHat = B*Z;

sHat(1,:) = sHat(1,:)/norm(sHat(1,:));
sHat(2,:) = sHat(2,:)/norm(sHat(2,:));


if norm(sHat(1,:)-s1) >= 0.5 && norm(-sHat(1,:)-s1) >= 0.5
sHat = [sHat(2,:);sHat(1,:)];
end
if norm(sHat(1,:)-s1) >= norm(-sHat(1,:)-s1)
    sHat(1,:) = -sHat(1,:);
end
if norm(sHat(2,:)-s2) >= norm(-sHat(2,:)-s2)
    sHat(2,:) = -sHat(2,:);
end
Error = sHat-[s1;s2];
normError = norm(Error);
ErrorPd = ErrorPd + (normError/norm([s1;s2]))^2;
end
figure(6)
ErrorPd = ErrorPd/100
plot(s1)
hold on 
plot(sHat(1,:))
legend('S1','sHat1')
title('Compare S1 With noisy sHat1 ')
figure(7)
plot(s2)
hold on 
plot(sHat(2,:))
legend('S2','sHat2')
title('Compare S2 With noisy sHat2 ')
%%%
B = [0 0 ;0 0] ;  
ErrorPe = [];
for K = 2:5
Esum = 0;
for itr = 1:100
 W = randn(size(X));
 W = W/norm(W);
 sig = sqrt(norm(X)^2/100);
 Y = X+sig*W;
 [U,D] = eig(Y*Y.');
 Z = sqrt(inv(D))*U.'*Y;
 for i = 1 : 2
   Ri = zeros(2);
  for k = 1:K
     ind = setdiff(1:2,i);
     Zi = Z(:,20*k - 19:20*k);
     Rzi = Zi*Zi.';
     Ri = Ri +(Rzi*B(ind,:).'* B(ind,:) *Rzi.');
  end  
     [V,lambdai] = eig(Ri,'vector');
     [lambdai,indi] = sort(lambdai);
     V = V(:,indi);   
     B(i,:) = V(:,1);
     
  if i > 1
     Rb = B(1:i-1,:).'*B(1:i-1,:);
      B(i,:) = (eye(2)-Rb)*(B(i,:).');
      B(i,:) = B(i,:)/norm(B(i,:));
  end
 end
sHat = B*Z;

sHat(1,:) = sHat(1,:)/norm(sHat(1,:));
sHat(2,:) = sHat(2,:)/norm(sHat(2,:));


if norm(sHat(1,:)-s1) >= 0.5 && norm(-sHat(1,:)-s1) >= 0.5
sHat = [sHat(2,:);sHat(1,:)];
end
if norm(sHat(1,:)-s1) >= norm(-sHat(1,:)-s1)
    sHat(1,:) = -sHat(1,:);
end
if norm(sHat(2,:)-s2) >= norm(-sHat(2,:)-s2)
    sHat(2,:) = -sHat(2,:);
end
Error = sHat-[s1;s2];
normError = norm(Error);
Esum = Esum + (normError/norm([s1;s2]))^2;
end
ErrorPe(k-1) = Esum/100;
end
ErrorPe
figure(8)
plot(2:5,ErrorPe)
title('decreasing Mean Error with increase windows number');
ylabel('Mean Error')
xlabel('k');



 %%%%%
W = randn(size(X));
W = W/norm(W);

[U,D] = eig(Y*Y.');
Z = sqrt(inv(D))*U.'*Y;
B = [0 0 ;0 0] ;  
ErrorPf = [];

for j = 1:4
  Esum = 0;
for itr = 1:100
W = randn(size(X));
W = W/norm(W);
sig = sqrt(norm(X)^2/(sqrt(10)^j));
Y = X+sig*W;
[U,D] = eig(Y*Y.');
Z = sqrt(inv(D))*U.'*Y;
 for i = 1 : 2
   Ri = zeros(2);
  for k = 1:K
     ind = setdiff(1:2,i);
     Zi = Z(:,20*k - 19:20*k);
     Rzi = Zi*Zi.';
     Ri = Ri +(Rzi*B(ind,:).'* B(ind,:) *Rzi.');
  end  
     [V,lambdai] = eig(Ri,'vector');
     [lambdai,indi] = sort(lambdai);
     V = V(:,indi);   
     B(i,:) = V(:,1);
     
  if i > 1
     Rb = B(1:i-1,:).'*B(1:i-1,:);
      B(i,:) = (eye(2)-Rb)*(B(i,:).');
      B(i,:) = B(i,:)/norm(B(i,:));
  end
 end
sHat = B*Z;
Hat(1,:) = sHat(1,:)/norm(sHat(1,:));
sHat(2,:) = sHat(2,:)/norm(sHat(2,:));


if norm(sHat(1,:)-s1) >= 0.5 && norm(-sHat(1,:)-s1) >= 0.5
sHat = [sHat(2,:);sHat(1,:)];
end
if norm(sHat(1,:)-s1) >= norm(-sHat(1,:)-s1)
    sHat(1,:) = -sHat(1,:);
end
if norm(sHat(2,:)-s2) >= norm(-sHat(2,:)-s2)
    sHat(2,:) = -sHat(2,:);
end
Error = sHat-[s1;s2];
normError = norm(Error);
Esum = Esum + (normError/norm([s1;s2]))^2;
 end
ErrorPf(j) = Esum/100;
end
ErrorPf
figure(9)
plot(5:5:20,ErrorPf)
title('decreasing Mean Error with increase SNR');
ylabel('Mean Error')
xlabel('SNR');
