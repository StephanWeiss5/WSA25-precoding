% Generate_Figures45.m
%
% Generate Figures 4 and 5 for the paper:
%   M.A. Bakhit, F.A. Khattak, S.J. Schlecht, G.W. Rice, and S. Weiss: 
%   "Challenges to Subcarrier MIMO Precoding and Equalisation with Smooth 
%   Phase Responses," 28th Workshop on Smart Antennas, Erlangen, Germany, 
%   September 2025.

clear all; close all;

U = zeros(2,2,2);
U(1,:,1) = [1,1];
U(2,:,2) = [1,-1];
U = U/sqrt(2);

Sigma = zeros(2,2,3);
%Sigma(1,1,:) = [.125*(1-1i) 1 .125*(1+1i)];
Sigma(1,1,:) = [1i 2 -1i]/2;
Sigma(2,2,:) = [0 1 1];

V = zeros(2,2,1);
V(:,:,1) = dftmtx(2)/sqrt(2);

A = PolyMatConv(U,PolyMatConv(Sigma,ParaHerm(V)));
% squared Euclidean norm is 1.75

%-----------------------------------------------------------------
% measuring distributions of singular values
%-----------------------------------------------------------------
randn('seed',0); Nfft = 2^12;
E = randn(size(A))+1i*randn(size(A))/sqrt(2)*1.75;
Ahat = A + E*0.01;
Ahatf = fft(Ahat,Nfft,3);
U1 = zeros(2,Nfft); U2 = zeros(2,Nfft); V1 = zeros(2,Nfft); V2 = zeros(2,Nfft);
SS = zeros(2,Nfft);
for k = 1:Nfft,
   [u,s,v] = svd(Ahatf(:,:,k));
   U1(:,k) = u(:,1); U2(:,k) = u(:,2);       
   V1(:,k) = v(:,1); V2(:,k) = v(:,2);       
   SS(:,k) = diag(s);
   if k == 1;
      angleU(1,k) = abs(U1(:,k)'*U1(:,k));
      angleU(2,k) = abs(U2(:,k)'*U1(:,k));
      angleV(1,k) = abs(V1(:,k)'*V1(:,k));
      angleV(2,k) = abs(V2(:,k)'*V1(:,k));
   else
      angleU(1,k) = abs(U1(:,k)'*U1(:,1));
      angleU(2,k) = abs(U2(:,k)'*U1(:,1));
      angleV(1,k) = abs(V1(:,k)'*V1(:,1));
      angleV(2,k) = abs(V2(:,k)'*V1(:,1));
   end;      
end;

%------------------------------------------------------------------------------
%  Figure 4
%------------------------------------------------------------------------------
figure(4);
f2 = (0:Nfft-1)/Nfft;
% dummies for legend
plot([-1 -1],[-1 -1],'b-'); hold on; plot([-1 -1],[-1 -1],'r--');
% now for serious
S1 = abs(fft([1i 2 -1i]/2,512)); S2 = abs(fft([0 1 1],512)); f = (0:511)/512;
h = plot(f,abs(S1),'-','color',[1 1 1]*0.8,'linewidth',3); 
h = plot(f,abs(S2),'-','color',[1 1 1]*0.8,'linewidth',3); 
plot(f2,SS(1,:),'b-');
plot(f2,SS(2,:),'r--');
plot([1 1]/2,[0 1.1],'k-.');
plot([1 1]/2,[0 1],'ko','MarkerFaceColor','k');
axis([0 1 0 2.1]); grid on;
legend({'$m=1$','$m=2$'},'interpreter','latex','location','NorthWest');
xlabel('norm.~angular freq.~$\Omega$','interpreter','latex');
ylabel('$\hat{\sigma}_m(\mathrm{e}^{\mathrm{j}\Omega})$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$',...
      '$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',(0:0.5:3),'YTickLabel',...
     {'$0$','$\frac12$','$1$','$\frac32$','$2$'});

%-----  insert #1
Box1 = [52/128 55/128 7/16 9/16];
plot(Box1([1 1 2 2 1]),Box1([3 4 4 3 3]),'k-','linewidth',1);
Box2 = [62/128 66/128 0 1/8];
plot(Box2([1 1 2 2 1]),Box2([3 4 4 3 3]),'k-','linewidth',1);
axes('position',[.4 .65 .1 .25]);
   h = plot(f,abs(S1),'-','color',[1 1 1]*0.8,'linewidth',3); 
   hold on;
   h = plot(f,abs(S2),'-','color',[1 1 1]*0.8,'linewidth',3); 
   plot(f2,SS(1,:),'b-');
   plot(f2,SS(2,:),'r--');
   axis(Box1);
    set(gca,'TickLabelInterpreter','latex',...
         'XTick',[52/128 55/128], 'XTickLabel',{'$\frac{13}{16}\pi$','$\frac{55}{64}\pi$'},...
         'YTick',[7 9]/16, 'YTickLabel',{'$\frac{7}{16}$','$\frac{9}{16}$'});
   grid on;
   
%-----  insert #2
axes('position',[.68 .25 .2 .2]);
   h = plot(f,abs(S1),'-','color',[1 1 1]*0.8,'linewidth',3); 
   hold on;
   h = plot(f,abs(S2),'-','color',[1 1 1]*0.8,'linewidth',3); 
   plot(f2,SS(1,:),'b-');
   plot(f2,SS(2,:),'r--');
   axis(Box2);
    set(gca,'TickLabelInterpreter','latex',...
         'XTick',[62 64 66]/128, 'XTickLabel',{'$\frac{31}{32}\pi$','$\pi$','$\frac{33}{32}\pi$'},...
         'YTick',[0 1]/8, 'YTickLabel',{'$0$','$\frac18$'});
   grid on;


%set(gcf,'OuterPosition',[230 250 570 300]);
set(gcf,'OuterPosition',[230 250 570 348]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc WSA25_2Fig4.eps

%------------------------------------------------------------------------------
%  Figure 5
%------------------------------------------------------------------------------
figure(5);
% dummies for legend
plot([-1 -1],[-1 -1],'b-'); hold on; plot([-1 -1],[-1 -1],'r--');
% for serious
plot((0:1024)/1024,abs(cos(pi*(0:1024)/1024)),'-','color',[1 1 1]*0.85,'linewidth',3); 
hold on;
plot((0:1024)/1024,abs(sin(pi*(0:1024)/1024)),'--','color',[1 1 1]*0.85,'linewidth',3); 
plot((0:Nfft-1)/Nfft,angleU(1,:),'b-');
plot((0:Nfft-1)/Nfft,angleU(2,:),'r--');
axis([0 1 0 1]); grid on;
legend({'$m=1$','$m=2$'},'interpreter','latex','location','SouthWest');
xlabel('norm.~angular freq.~$\Omega$','interpreter','latex');
ylabel('$\alpha_m(\Omega)$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$',...
      '$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',(0:0.25:1),'YTickLabel',...
     {'$0$,','$\frac{\pi}{8}$','$\frac{\pi}{4}$','$\frac{3\pi}{8}$','$\frac{\pi}{2}$'});
%set(gcf, 'OuterPosition', [230 250 570 320]);
set(gcf, 'OuterPosition', [230 250 570 370]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
print -depsc WSA25_2Fig5.eps


