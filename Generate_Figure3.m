% Generate_Figure3.m
%
% Generate Figure 3 for the paper:
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

%-----------------------------------------------------------------
% measuring distributions of singular values
%-----------------------------------------------------------------
W = ones(2,2,4);    % this is used to calculate the DFT coefficient at Omega = pi
W(:,:,[1 3]) = ones(2,2,2);
W(:,:,[2 4]) = -ones(2,2,2);
if exist('WSA25_2Fig3.mat')~=2,
   ss = zeros(2,10000);
   for i = 1:10000,
      E = randn(size(A))+1i*randn(size(A));
      Ahat = A + E*0.01;
      Ahatf = sum(Ahat.*W,3);        % Fourier coefficient at Omega = pi
      [~,s,~] = svd(Ahatf);
      ss(:,i) = diag(s);
      if mod(i,100)==0, disp(sprintf('%d iterations completed',i)); end;
   end;   
   save WSA25_2Fig3.mat ss
else
   load WSA25_2Fig3.mat   
end;
   
FS = 10;
set(0, 'DefaultTextInterpreter', 'latex', ...
       'DefaultAxesTickLabelInterpreter', 'latex', ...
       'DefaultLegendInterpreter', 'latex', ...
       'DefaultAxesFontSize', FS, ...
       'DefaultTextFontSize', FS);

figure(1);
if exist('WSA25_2Fig3a.mat')~=2,
   h = histfit(ss(1,:),100,'rician'); 
   Rice1X = h(2).XData;
   Rice1Y = h(2).YData;
   save WSA25_2Fig3a.mat Rice1X Rice1Y   
else
   load WSA25_2Fig3a.mat
   hist(ss(1,:),(0.94:0.0012:1.06)); hold on;
   plot(Rice1X,Rice1Y*0.68,'r-','linewidth',2);
end;   

xlabel('$\varsigma_1$','interpreter','latex'); ylabel('$\hat{p}(\varsigma_1)$','interpreter','latex');
axis([0.94 1.06 0 400]); grid on;
text(0.941,350,'(b)','interpreter','latex');
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0.94:0.02:1.06),'XTickLabel',{'$.94$','$.96$','$.98$','$1$','$1.02$','$1.04$','$1.06$'},...
    'YTick',(0:200:400),'YTickLabel',...
     {'$0$,','$200$','$400$'});
%set(gcf,'OuterPosition',[230 250 300 220]);
set(gcf,'OuterPosition',[230 250 300 285]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc WSA25_2Fig3a.eps

figure(2);
if exist('WSA25_2Fig3b.mat')~=2,
   h = histfit(ss(2,:),45,'rician');
   Rice2X = h(2).XData;
   Rice2Y = h(2).YData;
   save WSA25_2Fig3b.mat Rice2X Rice2Y   
else
   load WSA25_2Fig3b.mat
   hist(ss(2,:),(0:0.0012:0.0725)); hold on;
   plot(Rice2X,Rice2Y*0.68,'r-','linewidth',2);
end;   

xlabel('$\varsigma_2$','interpreter','latex'); ylabel('$\hat{p}(\varsigma_2)$','interpreter','latex');
axis([0 0.07 0 400]); grid on;
text(0.001,350,'(a)','interpreter','latex');
set(gca,'TickLabelInterpreter','latex',...
    'XTick',([0 0.02 0.04 0.06]),'XTickLabel',{'$0.0$','$0.02$','$0.04$','$.06$'},...
    'YTick',(0:200:600),'YTickLabel',...
     {'$0$,','$200$','$400$'});

%set(gcf,'OuterPosition',[230 250 195 220]);
set(gcf,'OuterPosition',[230 250 195 285]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc WSA25_2Fig3b.eps

