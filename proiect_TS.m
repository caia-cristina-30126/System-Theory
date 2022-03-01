
%% stabilitate INTERNA
A=[-1.38 -0.46 1.38;
    -0.92 -0.92 2.77;
    -0.92 -0.92 -0.92]
B=[1.38; 0.92; 0.92]; 
C=[-1 -1 -1]; 
D=[1]; 
sys=ss(A,B,C,D);
t=0:0.01:25 ;
u=(t>=0);
[y,t,x]=lsim(sys,u,t);
plot(t,x);

%% stabilitate EXTERNA
A=[-1.38 -0.46 1.38;
    -0.92 -0.92 2.77;
    -0.92 -0.92 -0.92]
B=[1.38; 0.92; 0.92]; 
C=[-1 -1 -1]; 
D=[1]; 
sys=ss(A,B,C,D);
t=0:0.01:25 ;
u=(t>=0);x0=[0 6 2];
[y,t,x]=lsim(sys,u,t,x0);
plot(t,y);

%% Lyapunov
A=[-1.38 -0.46 1.38;
    -0.92 -0.92 2.77;
    -0.92 -0.92 -0.92]
B=[1.38; 0.92; 0.92]; 
C=[-1 -1 -1]; 
D=[1];

Q=eye(3);
P=lyap(A',Q)
Val_pr_P=eig(P)

%%simulare pt conditii initiale
t=0:1:25;
	%u=(t>=0)
u=zeros(1,length(t));
x0=[11 8 1];
[y,t,x]=lsim(ss(A,B,C,D),u,t,x0);
figure;
V=zeros(1,length(t));
for i=1:length(t)
    V(i)=x(i,:)*P*x(i,:)';
end
plot(t,V);grid
eig(A);
Q=2*eye(3);
P =lyap(A',Q)
eig(P) 

%% raspunsul la impuls
H=tf([1 0 0 0],[1 3.23 6.84 3.17])
impulse(H)
hold on
t=0:0.1:10
h=dirac(t)-0.05*exp(-0.6*t)-3.18*exp(-1.31*t).*(cos(1.87*t)+0.1*sin(1.87*t))
plot(t,h,'r*')
%% raspunsul la treapta
H=tf([1 0 0 0],[1 3.23 6.78 3.17])
step(H)
hold on
t=0:0.1:10
y=0.08*exp(-0.6*t)+exp(-1.31*t).*(0.92*cos(1.87*t)-1.01*sin(1.87*t))
plot(t,y,'r*')
%% Raspunsul la rampa
H=tf([1 0 0 0],[1 3.23 6.84 3.17])
t=0:0.1:10
u=t %rampa
lsim(H,u,t);
hold on
y=-0.14*exp(-0.6*t)+0.14*exp(-1.31*t).*(cos(1.87*t)+3.95*sin(1.87*t))
plot(t,y,'r*')

%% parametrii Markov
num = [1 0 0 0];
den = [1 3.23 6.84 3.17];
markov = deconv([num zeros(1,6)],den)

%% FT minim
num = [1 0 0 0];
den = [1 3.23 6.78 3.17];


H=tf(num,den)
Hm=minreal(H)

%% FCC

A=[-3.23 -6.84 -3.17;
    1 0 0;
    0 1 0];
b=[1; 0; 0]; 
c=[-3.23 -6.84 -3.17]; 
d=[1]; 
num = [1 0 0 0];
den = [1 3.23 6.78 3.17];
b2=-3.23; b1=-6.48; b0=-3.17;
a2=3.23; a1=6.48; a0=3.17;
[A,b,c,d]=tf2ss(num,den); %obtinerea FCC (FT to ss)
   
sistem = ss(A,b,c,d);
x0=[0,0,0] % conditii initiale
t=0:0.01:5; 
u=ones(1,length(t)); % semnal treapta
[y,t,x] = lsim(sistem,u,t,x0); %raspunsul sistemului afisat grafic

subplot(121); plot(t,y); legend('y'); grid;
subplot(122); plot(t,x); legend('x.1', 'x.2', 'x.3'); grid;

%% FCO -

A=[-3.23 -6.84 -3.17;
    1 0 0;
    0 1 0];
b=[1; 0; 0]; 
c=[-3.23 -6.84 -3.17]; 
d=[1]; 
num = [1 0 0 0];
den = [1 3.23 6.78 3.17];
b2=-3.23; b1=-6.48; b0=-3.17;
a2=3.23; a1=6.48; a0=3.17;
[A,b,c,d]=tf2ss(num,den); %obtinerea FCC (FT to ss)
Afco = A'
Bfco=c'
Cfco=b'
Dfco=d
sistem = ss(Afco,Bfco,Cfco,Dfco);
x0=[9,0,2] % conditii initiale
t=0:0.01:5; 
u=ones(1,length(t)); % semnal treapta
[y,t,x] = lsim(sistem,u,t,x0); %raspunsul sistemului afisat grafic

subplot(121); plot(t,y); legend('y'); grid;
subplot(122); plot(t,x); legend('x.1', 'x.2', 'x.3'); grid;



%%
P=[1 3.23 6.78 3.12];
roots(P)

%% - eroarea stationara la pozitie
num = [1.62];
den = [1 2.61 5.25];
H=tf(num,den);
step(H);

%% 9. PERFORMANTELE
num = [1];
den = [1 3.23 6.84 3.17];
H=tf(num,den);
step(H,80);

%% trasarea LR
H=tf([1 0 0 0],[1 3.23 6.78 3.17])
rlocus(H)
rltool(H)

%% k apropiere
H=tf([1 0 0 0],[1 3.23 6.78 3.17])
-1/evalfr(H,0)

%% 11 c - trasarea

H=tf([5000 8000 9690 20520 9510], [3000 9690 20520 9510 0])
%rlocus(H)
% rltool(H)

%% 11 c

p1=[5000 8000 9690 20520 9510];
p2=[3000 9690 20520 9510 0];
roots(p1)
roots(p2)

%% unghiurile de plecare din poli
x1=atan2d(-4,3)+atan2d(4,3)-atan2d(0,0.000001)-2*atan2d(0,2.66)


%% unghiurile de sosire
x1=-atan2d(-8,0)+2*atan2d(-4,-3)+2*atan2d(-4,-0.34)


%% puncte de intalnire
p1=[0 20e3 16e2 2*9690 2052]
p2=[3e3 9690 20520 9510 0]
A1=conv(p1,p2)

conv([1 0],conv([1 5],conv([1 6], [1 2 2])))
p3=[5e3 8e3 9690 20520 9510]
p4=[0 12e3 3*9690 2*20520 9510]
A2=conv(p3,p4)

P=A1-A2
roots(P)



 


