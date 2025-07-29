% Generate_Figures12.m
%
% Generate Figures 1 and 2 for the paper:
%   M.A. Bakhit, F.A. Khattak, S.J. Schlecht, G.W. Rice, and S. Weiss: 
%   "Challenges to Subcarrier MIMO Precoding and Equalisation with Smooth 
%   Phase Responses," 28th Workshop on Smart Antennas, Erlangen, Germany, 
%   September 2025.

clear all; close all;

U = zeros(2,2,2);
U(1,:,1) = [1,1]/sqrt(2);
U(2,:,2) = [1,-1]/sqrt(2);
U = U/sqrt(2);

Sigma = zeros(2,2,3);
%Sigma(1,1,:) = [.125*(1-1i) 1 .125*(1+1i)];
Sigma(1,1,:) = [1i 2 -1i]/2;
Sigma(2,2,:) = [0 1 1];

V = zeros(2,2,1);
V(:,:,1) = dftmtx(2)/sqrt(2);

A = PolyMatConv(U,PolyMatConv(Sigma,ParaHerm(V)));

s1 = squeeze(Sigma(1,1,:)); s2 = squeeze(Sigma(2,2,:)); 
S1 = fft(s1,1024).*exp(1i*2*pi/1024*(0:1023)');
S2 = fft(s2,1024).*exp(1i*3*pi/1024*(0:1023)');


%----------------------------------------------------------
%  Figure 1
%----------------------------------------------------------
set_figure_style;
figure(1);
A2 = A*sqrt(2);
subplot(221);
stem(0:3,real(squeeze(A2(1,1,:))),'b'); hold on;
plot(0:3,imag(squeeze(A2(1,1,:))),'r*');
axis([-.5 3.5 -.6 1.2]); grid on;
ylabel('$\sqrt{2} C_{1,m}[n]$'); title('$m=1$');
subplot(222);
stem(0:3,real(squeeze(A2(1,2,:))),'b'); hold on;
plot(0:3,imag(squeeze(A2(1,2,:))),'r*');
axis([-.5 3.5 -.6 1.2]); grid on;
title('$m=2$');
subplot(223);
stem(0:3,real(squeeze(A2(2,1,:))),'b'); hold on;
plot(0:3,imag(squeeze(A2(2,1,:))),'r*');
axis([-.5 3.5 -.6 1.2]); grid on;
ylabel('$\sqrt{2} C_{2,m}[n]$'); xlabel('discrete time $n$');
subplot(224);
stem(0:3,real(squeeze(A2(2,2,:))),'b'); hold on;
plot(0:3,imag(squeeze(A2(2,2,:))),'r*');
axis([-.5 3.5 -.6 1.2]); grid on;
xlabel('discrete time $n$');
set(gcf, 'OuterPosition', [230 250 570 350]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
print -depsc WSA25_2Fig1.eps

%----------------------------------------------------------
%  Figure 2
%----------------------------------------------------------
figure(2);
Omega = (0:1024)/1024;
s1 = 1-sin(4*pi*Omega);
s2 = 2*cos(2*pi*Omega);
plot(Omega,s1,'b-'); hold on;
plot(Omega,s2,'r--');
plot(Omega(256:768),abs(s2(256:768)),'-.','color',[0 0.5 0.5]);
plot([0.25 0.75],[0 0],'ko');
axis([0 1 -2 2]); grid on;
set(gca, 'TickLabelInterpreter', 'latex', ...
        'XTick', (0:1/8:1), 'XTickLabel', {'$0$', '$\frac{\pi}{2}$', '$\pi$',...
        '$\frac{3\pi}{2}$', '$2\pi$', '$\frac{5\pi}{2}$', '$3\pi$',...
        '$\frac{7\pi}{2}$', '$4\pi$'}, ...
        'YTick', (-2:1:2), 'YTickLabel', {'$-2$', '$-1$', '$0$', '$1$', '$2$'});
legend({'$m=1$','$m=2$','$|\sigma_2^\prime(\Omega)|$'},'location','SouthEast');
ylabel('$\sigma_m^\prime(\Omega)$');
xlabel('norm.~angular freq. $\Omega$');
set(gcf, 'OuterPosition', [230 250 570 300]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
print -depsc WSA25_2Fig2.eps

return;

figure(5);
plot((0:1023)/1024,abs(S1));
hold on;
plot((0:1023)/1024,abs(S2),'r--');
S = abs(real([S1 S2].'));
S = sort(S,'descend');
S = [S; -S(end,:)];
Sdmin = min(abs(diff(S)));
figure(3);
subplot(211);
plot((0:1023)/1024,S');
subplot(212);
plot((0:1023)/1024,Sdmin');

