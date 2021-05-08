% LinearizedDynamicsOde.m
% 
% Description:
%   Evaluate the linearized dynamics (A, B)at the nominal state and control
%   vector.
%
% Inputs:
%   X - 4x1 state vector
%   u - 2x1 control vector
%   w - 4x1 process noise
% Outputs:
%   A_tilde - 4x4 linearized sys matrix, i.e. jacobian of nonlinear dynamics
%   B_tilde - 4x2 linearized ctrl matrix, i.e. jacobian of nonlinear dynamics

function [A_tilde, B_tilde] = LinearizedDynamicsOde(X, u, w) 
    
    global mu

%     A_11 = 0.0;
%     A_12 = 1.0;
%     A_13 = 0.0;
%     A_14 = 0.0;
% 
%     A_21 = 3*mu*(X(1)^2)/(((X(1)^2)+(X(3)^2))^(5/2)) -...
%         mu/(((X(1)^2)+(X(3)^2))^(3/2));
%     A_22 = 0.0;
%     A_23 = 3*mu*X(1)*X(3)/(((X(1)^2)+(X(3)^2))^(5/2));
%     A_24 = 0.0;
% 
%     A_31 = 0.0;
%     A_32 = 0.0;
%     A_33 = 0.0;
%     A_34 = 1.0;
% 
%     A_41 = 3*mu*X(1)*X(3)/(((X(1)^2)+(X(3)^2))^(5/2));
%     A_42 = 0.0;
%     A_43 =  3*mu*(X(3)^2)/(((X(1)^2)+(X(3)^2))^(5/2)) -...
%         mu/(((X(1)^2)+((X(3)^2)))^(3/2));
%     A_44 = 0.0;
% 
%     A_tilde = [A_11, A_12, A_13, A_14;
%             A_21, A_22, A_23, A_24;
%             A_31, A_32, A_33, A_34;
%             A_41, A_42, A_43, A_44];

    r = sqrt(X(1)^2 + X(3)^2);

    A_tilde = zeros(size(X,1), size(X,1));
    A_tilde(1,1) = 0.0;
    A_tilde(1,2) = 1.0;
    A_tilde(1,3) = 0.0;
    A_tilde(1,4) = 0.0;
    
    A_tilde(2,1) = 3.0*mu*(X(1)^2)/(r^5) - mu/(r^3);
    A_tilde(2,2) = 0.0;
    A_tilde(2,3) = 3.0*mu*X(1)*X(3)/(r^5);
    A_tilde(2,4) = 0.0;
    
    A_tilde(3,1) = 0.0;
    A_tilde(3,2) = 0.0;
    A_tilde(3,3) = 0.0;
    A_tilde(3,4) = 1.0;
    
    A_tilde(4,1) = 3.0*mu*X(1)*X(3)/(r^5);
    A_tilde(4,2) = 0.0;
    A_tilde(4,3) = 3.0*mu*(X(3)^2)/(r^5) - mu/(r^3);
    A_tilde(4,4) = 0.0;
    
    B_tilde = zeros(size(X,1), size(u,1));

    B_tilde(2,1) = 1.0;
    B_tilde(4,2) = 1.0;
    
    
end

    
    
    