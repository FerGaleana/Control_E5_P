function U = DOF_PLANAR_P2_u( X, t)

% Estados
x1 = X(1);  
x2 = X(2);  
x3 = X(3);  
x4 = X(4);

q1 = x1;
q2 = x2;
N = 5;
% Parametros del sistema
g = 9.81;
l1 = N; 
l2 = l1; 
lc1 = N/2; 
lc2 = lc1; 
m1 = N+1; 
m2 = m1;
b1 = 0.5; 
b2 = b1; 
I1 = 0.001*N; 
I2 = I1;

%Ecuaciones del sistema
d11 = (m1*(lc1^2))+m2*(l1^2+2*l1*lc2^2+2*l1*lc2*cos(q2))+I1+I2;
d12 = m2*((lc2^2)+l1*lc2*cos(q2))+I2;
d21 = d12;
d22 = m2*lc2^2+I2;
phi1 = (m1*lc1+m2*l1)*g*cos(q1+q2)+m2*lc2*g*cos(q1+q2);
phi2 = m2*lc2*cos(q1+q2);
c121 = -m2*l1*lc2*sin(q2);
c211 = c121;
c221 = c121;
c112 = -c121;

% Variable para guardar constante 
c = -m2*l1*lc2*sin(q2);
c1 = d22/(d11*d22-d12*d21);
c2 = d11/(d11*d22-d21*d12);
% Estados del sistema
f3 = c1*(-c121*x3*x4-c211*x4*x3-c221*(x4^2)-phi1+(d12/d22)*c112*(x3^2)+(d12/d22)*phi2);
f4 = c2*(-(d21/d11)*(-c121*x3*x4-c211*x4*x3-c221*(x4^2)-phi1)-c112*(x3^2)-phi2);
G = [c1,-c1*(d12/d22); -c2*(d21/d11),c2];
% Ganancias
Ke1 = -15;               % Ke1<0 polos estables
Ke2 = -15;               % Ke2<0 polos estables
Ke_ = [Ke1,Ke2];

K1 = Ke_(1);
K2 = Ke_(2);

Ref1 = (pi/3)*sin(N*t);
Ref2 = (pi/3)*cos(N*t);
dRef1 = (pi*N/3)*cos(N*t);   
dRef2 = -(pi*N/3)*sin(N*t);
ddRef1 = -(pi*(N^2)/3)*sin(N*t);   
ddRef2 = -(pi*(N^2)/3)*cos(N*t);

% Definición del error por bloques
% 1er bloque
e1 = [Ref1;  Ref2] - [x1; x2];
de1 = [dRef1; dRef2] - [x3; x4];
aux1 = [dRef1; dRef2] - K1*e1;
x3ref = aux1(1);
x4ref = aux1(2);

% 2do bloque
e2 = [x3ref; x4ref] - [x3; x4];
aux2 = [ddRef1; ddRef2] - K1*de1;
dx3ref = aux2(1);
dx4ref = aux2(2);

% Ley de control (bloques)
U = G^-1*([dx3ref; dx4ref] - [f3; f4] - K2*e2);
end