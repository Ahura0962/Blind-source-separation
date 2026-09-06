% be name khoda - Solution of Hw7
load('hw5.mat')
[M,N]=size(D);
N0=3;

%% a -- Subset Selection
Error=Inf(N,N,N);
for i=1:N-2
    for j=i+1:N-1
        for k=j+1:N
            Dsubset=D(:,[i j k]);
            shat=pinv(Dsubset)*x;
            Error(i,j,k)=norm(x-Dsubset*shat);
        end
    end
end

% Finding the minimum 
[val,pos]=min(Error(:));
khat=floor(pos/(N^2))+1;
pos1=pos-(khat-1)*N^2;
jhat=floor(pos1/(N))+1;
ihat=pos1-(jhat-1)*N;
possubset=[ihat jhat khat];
Dsubset=D(:,possubset);            %   [3     6     54] Selected Columns of D
s_subset=pinv(Dsubset)*x; %        %   [5     7     -3] Non-zeros values of s
% [val norm(x-Dsubset*s_subset)]

disp('Subset Section:')
[possubset; s_subset']
%% b --- norm 2
shat=pinv(D)*x;
stem(shat) % Check wheter shat is sparse or not!!!

%% c --- MP
x1=x;
posMP=zeros(1,N0);
sMP=zeros(N,1);
for i=1:N0
   ro=x1'*D;
   [~,posMP(i)]=max(abs(ro));
   sMP(posMP(i))=ro(posMP(i));
   x1=x1-sMP(posMP(i))*D(:,posMP(i));
end
disp('MP:')
[posMP; sMP(posMP)'] % Check wheter the solution is correct or not!!!

%% d --- OMP
x1=x;
posOMP=zeros(1,N0);
sOMP=zeros(N,1);
for i=1:N0
   ro=x1'*D;
   [~,posOMP(i)]=max(abs(ro));
   if i>1
       Dsub=D(:,posOMP(1:i));
       sOMP(posOMP(1:i))=pinv(Dsub)*x;
       x1=x-D*sOMP;
   else
       sOMP(posOMP(1))=ro(posOMP(1));
       x1=x1-sOMP(posOMP(1))*D(:,posOMP(1)); 
   end
end
disp('OMP:')
[posOMP; sOMP(posOMP)'] % Check wheter the solution is correct or not!!!

%% e --- BP
% Linear Programming
f=ones(2*N,1);
Aeq=[D -D];
beq=x;
lb=zeros(2*N,1);
yhat = linprog(f,[],[],Aeq,beq,lb,[]);
splus=yhat(1:N);
sminus=yhat(N+1:end);
sBP=splus-sminus;
posBP=find(abs(sBP)>0.01)';
disp('BP:')
[posBP;sBP(posBP)'] % Check wheter the solution is correct or not!!!

%% f --- IR
w=ones(N,1);
ITRmax=100;
for itr=1:ITRmax
    W=diag(w);
    y=pinv(D*(W^-1))*x;
    sIRLS=y./sqrt(w);
    for n=1:N
        if abs(sIRLS(n))<0.0001
            w(n)=1e10;
        elseif abs(sIRLS(n))>1e6
            w(n)=1e-10;
        else
            w(n)=1./abs(sIRLS(n));
        end
    end
end
posIRLS=find(abs(sIRLS)>0.1)';
disp('IRLS:')
[posIRLS;sIRLS(posIRLS)'] % Check wheter the solution is correct or not!!!

    