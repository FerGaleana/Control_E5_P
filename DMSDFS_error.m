function u = DMSDFS_error(X,t)
N = 5;
M1 = N; M2 = M1;
b1 = 3.5-0.3*N; b2 = b1;
k1 = 0.4*N; k2 = k1;
th = 10+7*N;
us1 = 0;   us2 = us1;
ud1 = 0; ud2 = ud1;
g = 9.8;

    A = [0            0           1             0; 
           0            0           0             1; 
        -(k1+k2)/M1    k2/M1    -(b1+b2)/M1     b2/M1;
        k2/M2          -k2/M2       b2/M2      -b2/M2];

    B = [0; 0; 0; 1/M2]; 
    C = [0.5*N 0 0 0];
    dx1 = A(1,:)*X;
    dx2 = A(2,:)*X;
    dx3 = A(3,:)*X - g*sin(th);
    dx4 = A(4,:)*X - g*sin(th);
    %Ganancias
    K1 = -20; K2 = -20; K3 = -20;
    %Ref. y deivadas
    yref = 4*sin(5*t) + 2.5*cos(t);
    dyref = 20*cos(5*t) - 2.5*sin(t);
    ddyref = -100*sin(5*t) - 2.5*cos(t);
    dddyref = -500*cos(5*t) + 2.5*sin(t);
    %Primer bloque
    e1 = yref - 2.5*X(1);
    de1 = dyref - 2.5*dx1;
    dde1 = ddyref - 2.5*dx3;
    x3ref = (1/2.5)*(dyref -K1*e1);
    dx3ref = (1/2.5)*(ddyref - K1*de1);
    ddx3ref = (1/2.5)*(dddyref - K1*dde1);
    %Segundo bloque
    e2 = x3ref - X(3);
    de2 = dx3ref - dx3;
    x4ref = (1/0.4)*(dx3ref + 0.8*X(1) - 0.4*X(2) + 0.8*X(3) - g*sin(th) - K2*e2);
    dx4ref = (1/0.4)*(ddx3ref + 0.8*dx1 - 0.4*dx2 + 0.8*dx3 - K2*de2);
   %Tercer bloque
   e3 = x4ref - X(4);
   u = (1/0.2)*(dx4ref - 0.4*X(1) + 0.4*X(2) - 0.4*X(3) + 0.4*X(4) - g*sin(th) - K3*e3);
      
end