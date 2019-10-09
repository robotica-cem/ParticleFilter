function whereIam(N)
   p=10*rand(2,N)-5;
   figure(1);plot(p(1,:),p(2,:),'r*');
   axis([-5 5 -5 5])
   figure(2);histogram2(p(1,:),p(2,:))
   pause
   
   p=0.5*randn(2,N)-1;
   figure(1);plot(p(1,:),p(2,:),'r*');
   axis([-5 5 -5 5])
   figure(2);histogram2(p(1,:),p(2,:))
   pause
   
   p=0.5*randn(2,round(3*N/4))+1;
   q=randn(2,round(N/4))-3;
   figure(1);plot(p(1,:),p(2,:),'r*',q(1,:),q(2,:),'r*');
   axis([-5 5 -5 5])
   figure(2);histogram2([p(1,:),q(1,:)],[p(2,:),q(2,:)])
end