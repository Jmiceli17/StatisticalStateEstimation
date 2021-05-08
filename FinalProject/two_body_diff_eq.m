% Nonlinear dynamics for satellite

function Xdot = two_body_diff_eq( t, X, u, w)
mu = 398600;
r = sqrt( X(1)^2 + X(3)^2  ); % x^2 + y^2
A = [ 0       1        0        0
    -mu/r^3 0        0        0
    0       0        0        1
    0       0        -mu/r^3  0 ];
B = [0 0; 1 0; 0 0; 0 1];
Gam = [0 0; 1 0; 0 0; 0 1];

if isempty(u) % No u, w
    
    u = [0 0]';
    
end
if isempty(w)
    
    w = [0 0]';
end
    
    Ap = [0,1,0,0
        3*mu*X(1)^2/((X(1)^2+X(3)^2)^(5/2)) - mu/((X(1)^2+X(3)^2)^(3/2)),0,3*mu*X(1)*X(3)/((X(1)^2+X(3)^2)^(5/2)),0
        0,0,0,1
        3*mu*X(1)*X(3)/((X(1)^2+X(3)^2)^(5/2)),0,3*mu*X(3)^2/((X(1)^2+X(3)^2)^(5/2)) - mu/((X(1)^2+X(3)^2)^(3/2)),0];
    Xdotp = Ap * X;
    
    % Standard nonlinear state space realization
    Xdot = A*X + Gam*w;
    
    if all(Xdotp ~= Xdot)
        keyboard;
    end
end