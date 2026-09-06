  
  X1w = X1(:,1:100);
  X1w2 = X1(:,34:133);
  Rx1 = X1w*X1w2.' ;
  X1w3 = X1(:,67:166);
  Rx2 = X1w*X1w3.' ;
  [Q,~] = eig(Rx1,Rx2);
  B = Q.';
  sHat = B*X1;
  fs = 100;
  t = 0 : length(sHat(1,:))-1 ;
  t = t/fs;
  figure(1)
  plot(t,sHat(1,:))
  title("Estimated source1")
  figure(2)
  plot(t,sHat(2,:))
  title("Estimated source2")
  figure(3)
  f = -fs/2 : fs/length(t): fs/2 -  fs/length(t);
  plot(f,abs(fftshift((fft(sHat(1,:))))));
  title("Furier Transform Estimated source1")
  figure(4)
  plot(f,abs(fftshift(fft(sHat(2,:)))));
  title("Furier Transform Estimated source2")
  %%%%% 
  X2w = X2(:,1:100);
  X2w2 = X2(:,34:133);
  Rx1 = X2w*X2w2.' ;
  X2w3 = X2(:,67:166);
  Rx2 = X2w*X2w3.' ;
  [Q,~] = eig(Rx1,Rx2);
  B = Q.';
  sHat = B*X2;
  figure(5)  
  plot(t,sHat(1,:))
  title("Estimated source1")
  figure(6)
  plot(t,sHat(2,:))
  title("Estimated source2")
  figure(7)
  plot(f,abs(fftshift(fft(sHat(1,:)))));
  title("Furier Transform Estimated source1")
  figure(8)
  plot(f,abs(fftshift(fft(sHat(2,:)))));
  title("Furier Transform Estimated source2")