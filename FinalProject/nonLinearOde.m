% nonLinearOde.m
% 
% Description:
%   Nonlinear system dynamics

function Xdot = nonLinearOde( t, X, u, w )
%     mu = 398600;
    global mu
    
    Gam = [0 0; 1 0; 0 0; 0 1];

    r = sqrt(X(1)^2 + X(3)^2);
    A = [ 0      1        0        0
        -mu/r^3  0        0        0
         0       0        0        1
         0       0        -mu/r^3  0 ];
%     Xdot = [X(2);
%             -mu*X(1)/(r^3);
%             X(4);
%             -mu*X(3)/(r^3)];

    Xdot = A*X + Gam*w;

end