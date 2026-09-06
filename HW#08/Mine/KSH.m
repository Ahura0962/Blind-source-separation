%% Part 1
clc
clear
close all
load("hw8.mat")
figure(1)
subplot(3,1,1);
plot(S(1,:))
title('1st source')
subplot(3,1,2);
plot(S(2,:))
title('2nd source')
subplot(3,1,3);
plot(S(3,:))
title('3th source')
X_Nonoise = A*S;
figure(2)
subplot(3,1,1);
plot(X_Nonoise(1,:))
title('1st observation')
subplot(3,1,2);
plot(X_Nonoise(2,:))
title('2nd observation ')
subplot(3,1,3);
plot(X_Nonoise(3,:))
title('3th observation ')
X = A*S + Noise;
figure(3)
subplot(3,1,1);
plot(X(1,:))
title('1st noisy observation')
subplot(3,1,2);
plot(X(2,:))
title('2nd noisy observation ')
subplot(3,1,3);
plot(X(3,:))
title('3th noisy observation ')
mu = 0.01;
B = rand(size(A));
B = normr(B);
error = inf;
T = 1001;
objective_func = zeros(1,6000);
for itr = 1:6000
    y = B*X;
    y_m = max(abs(y).').';
    y = y./y_m;
    B = B./y_m;
    Psi = zeros(3,T);
     for i = 1:size(y,1)
        ky = [ones(1,length(y(i,:)));y(i,:);y(i,:).^2;y(i,:).^3;y(i,:).^4;y(i,:).^5];
        EY = ky*ky.'/T;
        yPrime = [sum(zeros(1,size(y,2)));sum(ones(1,size(y,2)));sum(2*y(i,:));sum(3*y(i,:).^2);sum(4*y(i,:).^3);sum(5*y(i,:).^4)]/T;
        theta = (EY)^-1 *yPrime;
        Psi(i,:) = theta.'*ky;
     end
     dfdB = Psi*(X.')/T - ((B)^-1).';
     B = B - mu*dfdB;
     B = normr(B); 
     if norm(dfdB) < error
         B_optimum = B;
         error = norm(dfdB);
     end
     objective_func(itr) = norm(dfdB,'fro');
end
P = B_optimum*A;
Shat = B_optimum*X;

%% Part 2
[~ , sort_ind] = max(abs(P));
Pe = P(sort_ind,:);
Shat = Shat(sort_ind,:);
Shat = normr(Shat);
Shat(1,:) = Shat(1,:)/sign(Pe(1,1));
Shat(2,:) = Shat(2,:)/sign(Pe(2,2));
Shat(3,:) = Shat(3,:)/sign(Pe(3,3));
S = normr(S);
figure(4)
subplot(3,1,1);
plot(Shat(1,:))
hold on
plot(S(1,:))
legend('Estimated','Main')
title('1st source estimated with main source ')
subplot(3,1,2);
plot(Shat(2,:))
hold on
plot(S(2,:))
title('2nd source estimated with main source ')
subplot(3,1,3);
plot(Shat(3,:))
hold on 
plot(S(3,:))
title('3th source estimated with main source')
E = (norm(S-Shat ,'fro')/norm(S,'fro'))^2;
%% Part 3
figure(5)
plot(objective_func(1:2000))
title('objective func depends on itr')
