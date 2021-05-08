% LinearizedMeasurementOde.m
% 
% Description:
%   Evaluate the linearized sensing matrix at the nominal state and control
%   vector.
%
% Inputs:
%   X - 4x1 state vector
%   u - 2x1 control vector
%   v - 3x1 measurement noise vector
%   gsX - 4x1 state vector of the ground station
% Outputs:
%   C_tilde - 3x4 linearized measurement matrix, jacobian of nonlinear
%             measurments wrt X

function C_tilde = LinearizedMeasurementOde(X, u, v, gsX)

    % Initialize matrix
    C_tilde = zeros(3,4);



    rho = sqrt((X(1) - gsX(1))^2 + (X(3) - gsX(3))^2);  % range

%     C_tilde(1,1) = (X(1) - gsX(1))/rho;
%     C_tilde(1,2) = 0.0;
%     C_tilde(1,3) = (X(3) - gsX(3))/rho;
%     C_tilde(1,4) = 0.0;
% 
%         C_tilde(2,1) = (X(2) - gsX(2))/rho - ...
%             ( (X(1) - gsX(1))^2*(X(2)-gsX(2)) + (X(1) - gsX(1))*(X(3) - gsX(3))*(X(4) - gsX(4)) )/(rho^3);
% 
% %     C_tilde(2,1) = (rho*(X(2) - gsX(2)) - ((X(1) - gsX(1)) * (X(2) - gsX(2)) +...
% %         (X(3) - gsX(3))*(X(4) - gsX(4))) * C_tilde(1,1) ) / (rho^2);
% 
%     C_tilde(2,2) = (X(1) - gsX(1))/rho;
% 
%     C_tilde(2,3) = (X(4) - gsX(4))/rho - ...
%         ( (X(1) - gsX(1))*(X(2) - gsX(2))*(X(3)-gsX(3)) + (X(3) - gsX(3))^2*(X(4) - gsX(4)) )/(rho^3);
% 
% %         C_tilde(2,3) = (rho * (X(3) - gsX(3)) - ((X(1)-gsX(1)) *(X(2) - gsX(2)) +...
% %           (X(3) - gsX(3))*(X(4) - gsX(4))) * C_tilde(1,3)) / (rho^2);
% 
%     C_tilde(2,4) = (X(3) - gsX(3))/rho;
% 
%     C_tilde(3,1) = (gsX(3)-X(3))/( (gsX(1) - X(1))^2 + (gsX(3) - X(3))^2 );
%     C_tilde(3,2) = 0.0;
%     C_tilde(3,3) = (X(1) - gsX(1))/( (gsX(1) - X(1))^2 + (gsX(3) - X(3))^2 );
%     C_tilde(3,4) = 0.0;
    
    
    
    
    C_tilde(1,1) = (X(1) - gsX(1))/rho;
    C_tilde(1,2) = 0.0;
    C_tilde(1,3) = (X(3) - gsX(3))/rho;
    C_tilde(1,4) = 0.0;
    
    C_tilde(2,1) = (X(2) - gsX(2))/rho - ...
        ( (X(1) - gsX(1))^2*(X(2)-gsX(2)) + (X(1) - gsX(1))*(X(3) - gsX(3))*(X(4) - gsX(4)) )/rho^3;
    C_tilde(2,2) = (X(1) - gsX(1))/rho;
    C_tilde(2,3) = (X(4) - gsX(4))/rho - ...
        ( (X(1) - gsX(1))*(X(2) - gsX(2))*(X(3)-gsX(3)) + (X(3) - gsX(3))^2*(X(4) - gsX(4)) )/rho^3;
    C_tilde(2,4) = (X(3) - gsX(3))/rho;
    
    C_tilde(3,1) = -(X(3) - gsX(3))/( (gsX(1) - X(1))^2 + (gsX(3) - X(3))^2 );
    C_tilde(3,2) = 0.0;
    C_tilde(3,3) = (X(1) - gsX(1))/( (gsX(1) - X(1))^2 + (gsX(3) - X(3))^2 );
    C_tilde(3,4) = 0.0;
    
end

