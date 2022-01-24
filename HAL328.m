clear all;
close all;
clc;

%zadání:
a3 = 1;
a2 = 18.2;
a1 = 90.91;
a0 = 73.71;
b0 = 6.3;
citatel = [b0];
jmenovatel = [a3 a2 a1 a0];

% ad/ A
[A,B,C,D]=tf2ss(citatel,jmenovatel);

%pøevod na frobeniùv tvar:
FROA = rot90(A,2)
FROB = flipud(B)
FROC = fliplr(C)
FROD = D

% pøevod na Jordanùv tvar
[r,p,k]=residue(citatel,jmenovatel)
JORA=diag(p)
JORB=ones(size(p,1),1)
JORC=r'
JORD=0

figure
step(JORA,JORB,JORC,JORD)
title('h(t) - Jordanuv tvar')
figure
step(FROA,FROB,FROC,FROD)
title('h(t) - Frobeniuv tvar')

% ad/ B
syms s;
I = eye(3);
Fi_s = inv(s*I-JORA); % obraz matice prechodu
x_0 = [1;1.5;1]; % zvoleny poc. stav
x_t = ilaplace(Fi_s) * x_0
%vykresleni
figure   
hold on
ezplot(x_t(1),[0 6])
ezplot(x_t(2),[0 6])
ezplot(x_t(3),[0 6])
axis auto
legend('x_1(t)','x_2(t)','x_3(t)')  
title('Prubeh stavovych velicin')
grid on

% uplne reseni (vcetne buzeni)
U_s = 5/s  % obraz budiciho signalu - konstanta: 5*skok
X_0 = x_0; % ?? ./ s;
X_s_upl = Fi_s * (X_0 + JORB * U_s)
x_t_upl = ilaplace(X_s_upl)
%vykresleni
figure   
hold on
ezplot(x_t_upl(1),[0 6])
ezplot(x_t_upl(2),[0 6])
ezplot(x_t_upl(3),[0 6])
axis auto
legend('x_1(t)','x_2(t)','x_3(t)')
title('Prubeh stavovych velicin - buzena soustava')
grid on

% ad/ C
% riditelnost dosazitelnost systemu:
Q_rid = ctrb(JORA, JORB) % matice riditelnosti
rank_Q_rid = rank(Q_rid) % system je dosazitelny a riditelny, pokud rad matice odpovida radu systemu
% pozorovatelnost
Q_poz = obsv(JORA, JORC)
rank_Q_poz = rank(Q_poz)

% ad/ D
% reg
Gs = tf(citatel,jmenovatel);  % system
Gr = pidtune(Gs,'PID')   % navrhni PID reg - prenos
G0 = series(Gs,Gr); % prenos otevrene smycky
figure
margin(G0)

%Nyquistovo kritérium stability
figure
nyquist(G0) % na zaklade otevrene smycky

Gw = feedback(G0,1);% ?? ,-1); % uzavrit f.b.
stepinfo(Gw) % prekmit, doba nabehu ...
figure
step(Gw)
grid on

%Michaljovo kritérium stability
[numGw,denGw]=tfdata(Gw);
w=0:0.01:2000; % pozor uprava !!
Gjw=1./polyval(denGw{1,1},1i*w);
figure
plot(real(Gjw),imag(Gjw))
grid on;
ylabel('Im'), xlabel('Re')
title('Michailovo kriterium stability');

% ad/ E
Ge = feedback(1,G0,-1); % pro regulacni odchylku
[e_t,t] = step(Ge,10); % pro t = 0..20
dt = t(2)-t(1);   % casovy krok
% kvalita regulace
I_kvad = sum((e_t).^2.*dt) % int. kvadr. kriterium
ITAE = sum(abs(e_t).*t.*dt) % ITAE

% ad/ F
%uložení prùbìhu frekvenèní charakteristiky:
[ReGs,ImGs,omegGs] = nyquist(Gs,{0.01, 2000});
omegGs = omegGs';
ReGs1 = permute(ReGs, [1 3 2]);
ImGs1 = permute(ImGs, [1 3 2]);
% prechodova char z hodnot freq. - inv. Fourier
dt = 0.001;
time = 0:dt:8;
g_recon = zeros(length(time)); % alloc memory - faster
h_recon = zeros(length(time)); % alloc memory - faster
acc = 0; % accumulator
for k=1:length(time)
   value = 1/(2*pi)* (trapz(omegGs,(ReGs1 + i*ImGs1) .* exp(i*time(k)*omegGs)) + trapz(omegGs,(ReGs1 - i*ImGs1) .* exp(-i*time(k)*omegGs)));
   acc = acc + value * dt;
   g_recon(k) = value;
   h_recon(k) = acc;
end
figure
plot(time, h_recon,'g');
hold on
%plot(time, g_recon,'r');
step(Gs,time)

% ad/ G
Td = Gr.Kd/Gr.Kp
Ti = Gr.Kp/Gr.Ki
Tz = min( abs(roots(jmenovatel)))/10 % sampling ten-times faster than fastest time-constant
q0 = Gr.Kp*(1+Tz/Ti+Td/Tz)
q1 = -Gr.Kp*(1+2*Td/Tz)
q2 = Gr.Kp*Td/Tz
