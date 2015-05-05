% A final real-time ADMM 
% real-time applies to n16
%ADMM with attacks

clear all; clc; close all;
%load n16 % PST model  
%load IEEE68_dyn

%eigen valus conjugates
o_i=[... % eigenvalues
-0.3+2.2*1i,...
-0.35+3.2*1i,...
-0.4+3.5*1i];

%why residues
% residues
r=[-0.3+2.2*1i,-0.3+2.2*1i,-0.1+2.2*1i,...
-0.8+3.2*1i,-0.3+2.2*1i,-0.3+0.1*1i,...
-0.01+3.5*1i,-0.3+2.2*1i,-0.3+6.2*1i];

%initial A
actual_a=[4.6119,-9.6758,11.6141,-8.4050,3.4817,-0.6570];

T=0.2;
t=0:T:20;
y=zeros(10,length(t));

%theta
for i=1:10
    for j=1:3
        y(i,:)=y(i,:)+r(j)*exp(o_i(j)*t)+conj(r(j))*exp(conj(o_i(j)*t));
    end
end
        
Y1=y(1,:);        
Y2=y(2,:);
Y3=y(3,:);
Y4=y(4,:);
Y5=y(5,:);
Y6=y(6,:);
Y7=y(7,:);
Y8=y(8,:);
Y9=y(9,:);
Y10=y(10,:);


deg=6;

%hankel
m=40; % height of the Hankel matrix
H1=fliplr(Hankel(Y1,m,deg));
H2=fliplr(Hankel(Y2,m,deg));
H3=fliplr(Hankel(Y3,m,deg));
H4=fliplr(Hankel(Y4,m,deg));
H5=fliplr(Hankel(Y5,m,deg));
H6=fliplr(Hankel(Y6,m,deg));
H7=fliplr(Hankel(Y7,m,deg));
H8=fliplr(Hankel(Y8,m,deg));
H9=fliplr(Hankel(Y9,m,deg));
H10=fliplr(Hankel(Y10,m,deg));


% H1*x=C1, etc.
C1=Y1(deg+1:deg+m).';
C2=Y2(deg+1:deg+m).';
C3=Y3(deg+1:deg+m).';
C4=Y4(deg+1:deg+m).';
C5=Y5(deg+1:deg+m).';
C6=Y6(deg+1:deg+m).';
C7=Y7(deg+1:deg+m).';
C8=Y8(deg+1:deg+m).';
C9=Y9(deg+1:deg+m).';
C10=Y10(deg+1:deg+m).';

% 
w1=zeros(deg,1); w2=w1; w3=w1;w4=w1;w5=w1;w6=w1;w7=w1;w8=w1;w9=w1;w10=w1;
% 
% 
x1(:,1)=1*ones(deg,1);x2(:,1)=ones(deg,1);x3(:,1)=1*ones(deg,1); x4(:,1)=ones(deg,1);x5(:,1)=ones(deg,1);
z=transpose(mean([x1,x2,x3,x4,x5].'));
r=log(roots([1,-z(:,1).']));
% 
rho=10^-6;
N=40;
disp('-----------------------------------------------');
tic;
l=zeros(N,5);

%round robbin

%round robin strategy  -- g is the choice matrix

%10 estimators
%g1=[1;0;0;0;1;0;0;0;1];
gtemp=eye(10);
g1=gtemp(:);
g=[g1;g1;g1;g1;g1;g1;g1;g1;g1;g1;g1;g1];


TT=size(H1,1);
% Start of the iterations
flag=1;
for j=1:N
     
   %equation 10 ??
    x1(:,j+1)=(H1(1:TT,:)'*H1(1:TT,:)+rho*eye(size(H1(1:TT,:),2)))\(H1(1:TT,:)'*C1(1:TT,:)-w1+rho*z(:,j));
    x2(:,j+1)=(H2(1:TT,:)'*H2(1:TT,:)+rho*eye(size(H2(1:TT,:),2)))\(H2(1:TT,:)'*C2(1:TT,:)-w2+rho*z(:,j));
    x3(:,j+1)=(H3(1:TT,:)'*H3(1:TT,:)+rho*eye(size(H3(1:TT,:),2)))\(H3(1:TT,:)'*C3(1:TT,:)-w3+rho*z(:,j));
    x4(:,j+1)=(H4(1:TT,:)'*H4(1:TT,:)+rho*eye(size(H4(1:TT,:),2)))\(H4(1:TT,:)'*C4(1:TT,:)-w4+rho*z(:,j));
    x5(:,j+1)=(H5(1:TT,:)'*H5(1:TT,:)+rho*eye(size(H5(1:TT,:),2)))\(H5(1:TT,:)'*C5(1:TT,:)-w5+rho*z(:,j));
    x6(:,j+1)=(H6(1:TT,:)'*H6(1:TT,:)+rho*eye(size(H6(1:TT,:),2)))\(H6(1:TT,:)'*C6(1:TT,:)-w6+rho*z(:,j));
    x7(:,j+1)=(H7(1:TT,:)'*H7(1:TT,:)+rho*eye(size(H7(1:TT,:),2)))\(H7(1:TT,:)'*C7(1:TT,:)-w7+rho*z(:,j));
    x8(:,j+1)=(H8(1:TT,:)'*H8(1:TT,:)+rho*eye(size(H8(1:TT,:),2)))\(H8(1:TT,:)'*C8(1:TT,:)-w8+rho*z(:,j));
    x9(:,j+1)=(H9(1:TT,:)'*H9(1:TT,:)+rho*eye(size(H9(1:TT,:),2)))\(H9(1:TT,:)'*C9(1:TT,:)-w9+rho*z(:,j));
    x10(:,j+1)=(H10(1:TT,:)'*H10(1:TT,:)+rho*eye(size(H10(1:TT,:),2)))\(H10(1:TT,:)'*C10(1:TT,:)-w10+rho*z(:,j));

   %x1(:,j+1) = x1(:,j+1)+10;
%   x2(:,j+1) = x2(:,j+1)+5;
%  if j==7
%         x2(:,j+1) = x2(:,j+1)+4;
%  
%  end
  
  if flag==1
       flag=0;
   else if flag==0;
           flag=1;
       end
   end
   
    if j==1
        p=1;
    else 
        p=p+10;
    end
   
%    z(:,j+1)=transpose(mean([g(p)*x1(:,j+1),g(p+1)*x2(:,j+1),g(p+2)*x3(:,j+1),g(p+3)*x4(:,j+1),g(p+4)*x5(:,j+1),g(p+5)*x6(:,j+1),g(p+6)*x7(:,j+1),g(p+7)*x8(:,j+1),g(p+8)*x9(:,j+1),g(p+9)*x10(:,j+1)].'));
    z(:,j+1)=transpose(10*mean([g(p)*x1(:,j+1),g(p+1)*x2(:,j+1),g(p+2)*x3(:,j+1),g(p+3)*x4(:,j+1),g(p+4)*x5(:,j+1),g(p+5)*x6(:,j+1),g(p+6)*x7(:,j+1),g(p+7)*x8(:,j+1),g(p+8)*x9(:,j+1),g(p+9)*x10(:,j+1)].'));
    
     
%    z(:,j+1)=transpose(mean([x1(:,j+1),x2(:,j+1),x3(:,j+1),x4(:,j+1),x5(:,j+1),x6(:,j+1),x7(:,j+1),x8(:,j+1),x9(:,j+1),x10(:,j+1)].'));
    
%     l(j,1)=norm(H1(1:TT,:)*z(:,j+1)-C1(1:TT,:));
%     l(j,2)=norm(H2(1:TT,:)*z(:,j+1)-C2(1:TT,:));
%     l(j,3)=norm(H3(1:TT,:)*z(:,j+1)-C3(1:TT,:));
%     
%     l(j,:)=0.5*(l(j,:).*l(j,:));
    
    %step (e)
    w1=w1+rho*(x1(:,j+1)-z(:,j+1)); 
    w2=w2+rho*(x2(:,j+1)-z(:,j+1));
    w3=w3+rho*(x3(:,j+1)-z(:,j+1));
    w4=w4+rho*(x4(:,j+1)-z(:,j+1)); 
    w5=w5+rho*(x5(:,j+1)-z(:,j+1));
    w6=w6+rho*(x6(:,j+1)-z(:,j+1));
    w7=w7+rho*(x7(:,j+1)-z(:,j+1)); 
    w8=w8+rho*(x8(:,j+1)-z(:,j+1));
    w9=w9+rho*(x9(:,j+1)-z(:,j+1));
    w10=w10+rho*(x10(:,j+1)-z(:,j+1));
    
    
    r(:,j+1)=sort(log(roots([1,-z(:,j+1).']))/T);%convert them to root continuous time
    fprintf('%10.4f',log10(max(l(j,:))))

    
    %error
    for k=1:length(o_i)
        [i1,i2]=min(abs(o_i(k)-r(:,j+1)));
        m_real(k,j)=real(r(i2,j+1));
        m_imag(k,j)=imag(r(i2,j+1));
    end
    fprintf('\n');         

end


disp('-----------------------------------------------');
r1=log(roots([1;-pinv(H1)*C1]))/T;
r2=log(roots([1;-pinv(H2)*C2]))/T;
r3=log(roots([1;-pinv(H3)*C3]))/T;


rc=log(roots([1;-pinv([H1;H2;H3])*[C1;C2;C3]]))/T;
t=(0:T:T*(length(Y1)-1));


% h1=figure('Position', [100, 100, 800, 600]);
% plot(0:N-1,log10(l),'linewidth',2)
% set(gca,'FontSize',24)
% xlabel('Iteration (k)','Interpreter','latex','fontsize',24);
% ylabel('$\log_{10}(e_i)$','Interpreter','latex','fontsize',24);
% legend({'PDC 1','PCD 2','PDC 3'},'fontsize',20,'FontName','Times New Roman');
% xlim([0,N-1])
%saveas(gcf,'C:\Users\snabavi\Dropbox\TSG\figs2\pq1.eps','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h2=figure('Position', [100, 100, 800, 600]);
set(gca,'FontSize',24)
plot(0:N-1,z(:,2:end),'linewidth',2)
hold on
plot(0:N-1,repmat(actual_a,N,1),'--','Linewidth',2)
xlabel('Iteration (k)','Interpreter','latex','fontsize',24);
ylabel('$a$','Interpreter','latex','fontsize',24);
legend({'$a_1$','$a_2$','$a_3$','$a_4$','$a_5$','$a_6$'},'Interpreter','latex','fontsize',24,'Position',[0.75,0.63,0.14,0.28]);
xlim([0,N-1])
%saveas(gcf,'C:\Users\snabavi\Dropbox\TSG\figs2\pq2.eps','epsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h3=figure('Position', [100, 100, 800, 600]);
% plot(0:N-1,m_imag,'linewidth',2)
% hold on
% plot(0:N-1,repmat(imag(o_i),N,1),'--','Linewidth',2)
% set(gca,'FontSize',24)
% xlabel('Iteration (k)','Interpreter','latex','fontsize',24);
% ylabel({'$\Omega$ (rad/sec)'},'Interpreter','latex','fontsize',24);
% legend({'$\Omega_1$','$\Omega_2$','$\Omega_3$'},'Interpreter','latex','fontsize',24,'Position',[0.75,0.63,0.14,0.28]);
% xlim([0,N-1])
% ylim([2,10])
% %saveas(gcf,'C:\Users\snabavi\Dropbox\TSG\figs2\pq3.eps','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1=l;
% savefile='admm';
% save(savefile,'N','L1');

