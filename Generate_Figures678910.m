% Generate_Figures678910.m
%
% Generate Figures 6 to 10 for the paper:
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

randn('seed',0); Nfft = 2^12;
E = randn(size(A))+1i*randn(size(A))/sqrt(2)*1.75;
Ahat = A + E*0.001;

S = DinosaurBoneSVD_v4(Ahat);

%-------------------------------------------------------------------------------
%  Figure 6 minimum distance
%-------------------------------------------------------------------------------
load DinosaurBoneInterimFig6.mat
figure(6);
LD = length(D);
plot((0:LD-1)/LD,D,'b-');
hold on; plot([0 1],[1 1]*Thresh,'r--');
axis([0 1 0 1.5]); grid on;
xlabel('norm.~angular freq.~$\Omega$','interpreter','latex'); 
ylabel('$d_{\min}(\Omega)$, $\mathcal{T}$','interpreter','latex');
legend({'$d_{\min}(\Omega)$','$\mathcal{T}$'},'interpreter','latex','location','NorthEast');
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$',...
      '$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',(0:.5:1.5),'YTickLabel',...
     {'$0$','$0.5$','$1$','$1.5$'});
%set(gcf,'OuterPosition',[230 250 570 250]);
set(gcf,'OuterPosition',[230 250 570 280]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc WSA25_2Fig6.eps

%-------------------------------------------------------------------------------
%  Figure 7 segmentation
%-------------------------------------------------------------------------------
figure(7);
load DinosaurBoneInterimResults.mat
% plot ground truth
plot((0:1024)/1024, 1- sin(2*pi*(0:1024)/1024) , '-','linewidth',3,'color',[1 1 1]*0.85); hold on;
plot((0:1024)/1024, abs(2*cos(pi*(0:1024)/1024)) , '-','linewidth',3,'color',[1 1 1]*0.85);
% plot segments
for q = 1:3,
   plot( (Omega(q,1):Omega(q,2))/1024 , Sf(:,Omega(q,1):Omega(q,2))','b-');
   plot([Omega(q,1), Omega(q,1)]/1024,Sf(:,Omega(q,1))' ,'bo','MarkerFaceColor','b');
   plot([Omega(q,2), Omega(q,2)]/1024,Sf(:,Omega(q,2))' ,'bo','MarkerFaceColor','b');
end;
plot( (Omega(4,1):1024)/1024 , Sf(:,Omega(4,1):1024)' , 'b-');
plot( (0:Omega(4,2)-1024)/1024,Sf(:,1:Omega(4,2)-1023)','b-');
plot( [Omega(4,1), Omega(4,1)]/1024 , Sf(:,Omega(4,1))' ,'bo','MarkerFaceColor','b');
plot( [Omega(4,2)-1024, Omega(4,2)-1024]/1024 , Sf(:,Omega(4,2)-1023)' ,'bo','MarkerFaceColor','b');
axis([0 1 0 2]); grid on;
xlabel('norm.~angular freq.~$\Omega$','interpreter','latex'); 
ylabel('$\hat{\sigma}_m(\mathrm{e}^{\mathrm{j}\Omega})$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$',...
      '$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',(0:.5:2),'YTickLabel',...
     {'$0$','$0.5$','$1$','$1.5$','$2$'});
%set(gcf,'OuterPosition',[230 250 570 250]);
set(gcf,'OuterPosition',[230 250 570 280]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc WSA25_2Fig7.eps
   
%-------------------------------------------------------------------------------
%  Figure 8 partial reconstructions and alignment
%-------------------------------------------------------------------------------
figure(8);
%  M = 3; Q = 4;
M = size(St,1);
t = (-2:2);
Q = size(Omega,1);
CCode = [0 0 1; 1 0 0; 0 .5 0];
seq = [1 2 2 1   2 1 1 2]; 
% seq = mod(0:(M * Q - 1), size(CCode, 1)) + 1; % Ensure seq matches the subplot count
i = 0;
for m = 1:M,
   for q = 1:Q,
      i = i + 1;
      subplot(M,Q,(m-1)*Q+q); 
      h = plot([-2.45 -2.45 2.45 2.45 -2.45],[1.17 -1.17 -1.17 1.17 1.17],'-','linewidth',1);
      set(h(1),'color',CCode(seq(i),:)); hold on;
      stem(t,real(squeeze(St(m,:,q))),'b','linewidth',1); hold on; 
      stem(t,imag(squeeze(St(m,:,q))),'r*','linewidth',1);
      axis([-2.5 2.5 -1.2 1.2]); grid on;
   end;
end;

%%%%for legend only%
for i = 1:M,
    subplot(M,Q,(i-1)*Q+1);
    ylabel(sprintf('$s_{q,%d}[\\tau]$', i), 'Interpreter', 'latex','fontsize',18);
end 
for i = 1:Q,
    subplot(M,Q,i);
    title(sprintf('$q=$%d',i),'interpreter','latex','fontsize',16); 
end 
for i = 1:Q,
    subplot(M,Q,Q+i);
    xlabel('$\tau$','interpreter','latex','fontsize',18); 
end 
%set(gcf,'OuterPosition',[230 250 2*570 2*200]);
set(gcf,'OuterPosition',[230 250 2*570 2*250]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc WSA25_2Fig8.eps

%-------------------------------------------------------------------------------
%  Figure 9   reconstruction
%-------------------------------------------------------------------------------
figure(9);
% plot ground truth
plot((0:1024)/1024, 1- sin(2*pi*(0:1024)/1024) , '-','linewidth',3,'color',[1 1 1]*0.85); hold on;
% plot reconstruction
SS = zeros(2,1024); SS(:,1:3) = S(:,3:5); SS(:,1023:1024) = S(:,1:2);
Sf = fft(SS,1024,2);
f = (0:1023)/1024;
plot(f,real(Sf(2,:)),'b-');
plot(f,real(Sf(1,:)),'r--');
plot(f,imag(Sf(1,:)),'-.','color',[0 0.5 0.5]);
plot((0:1024)/1024, 2*cos(pi*(0:1024)/1024) , '-','linewidth',3,'color',[1 1 1]*0.85);
plot(f,real(Sf(1,:)),'r--');
axis([0 1 -2 2]); grid on;
xlabel('norm.~angular freq.~$\Omega$','interpreter','latex'); 
ylabel('singular values','interpreter','latex');
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$',...
      '$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',(-2:1:2),'YTickLabel',...
     {'$-2$','$-1$','$0$','$1$','$2$'});
legend({'ground truth','$s_1(\mathrm{e}^{\mathrm{j}\Omega})$',...
'$\Re\{s_2(\mathrm{e}^{\mathrm{j}\Omega})\}$','$\Im\{s_2(\mathrm{e}^{\mathrm{j}\Omega})\}$'},...
       'interpreter','latex','location','SouthEast');

%set(gcf,'OuterPosition',[230 250 570 250]);
set(gcf,'OuterPosition',[230 250 570 320]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc WSA25_2Fig9.eps

%-------------------------------------------------------------------------------
%  Figure 10   reconstruction of perturbed singular values
%-------------------------------------------------------------------------------
figure(10);
S = [ones(2,1)*eps, S, ones(2,1)*eps];
plot(-3:3,20*log10(abs(S(1,:))),'bo-');
hold on;
plot(-3:3,20*log10(abs(S(1,:))),'r*--');
text(-4.7,-25,'(a)','interpreter','latex');
axis([-5 5 -150 0]); grid on;
xlabel('time index $n$','interpreter','latex'); 
ylabel('$20\log_{10}\{|\hat{s}_m[n]|\}$','interpreter','latex');
%set(gcf,'OuterPosition',[230 250 150 250]);
set(gcf,'OuterPosition',[230 250 150 320]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc WSA25_2Fig10a.eps

figure(11);
Ahat = A + E*0.001;
Nfft = 2^15;
Ahatf = fft(Ahat,Nfft,3);
Shatf = zeros(2,Nfft);
for k = 1:Nfft,
   [~,s,~] = svd(squeeze(Ahatf(:,:,k)));
   Shatf(:,k) = diag(s);
end;
shat = fftshift(ifft(Shatf,Nfft,2).');
plot(-2^14:2^14-1,20*log10(abs(shat(:,1))),'b-');
hold on;
plot(-2^14:2^14-1,20*log10(abs(shat(:,2))),'r--');
text(-980,-25,'(b)','interpreter','latex');
axis([-1000 1000 -150 0]); grid on;
legend({'$m=1$','$m=2$'},'interpreter','latex','location','NorthEast');
xlabel('time index $n$','interpreter','latex'); 
ylabel('$20\log_{10}\{|\hat{\sigma}_m[n]|\}$','interpreter','latex');
%set(gcf,'OuterPosition',[230 250 420 250]);
set(gcf,'OuterPosition',[230 250 420 320]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc WSA25_2Fig10b.eps

