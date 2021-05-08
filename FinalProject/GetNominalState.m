% Get the nominal state of the satellite at this time
% (equivalent to evaluating symbolic xnom)


function x_nom = GetNominalState(time)

    global r0 n vel0

    x_nom = [r0*cos(n*time); 
            -vel0*sin(n*time); 
            r0*sin(n*time); 
            vel0*cos(n*time)];



end