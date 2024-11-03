

%% write a code in matlab to solve lamberts problem
function [V1,V2]=lambert(mu,R0,R1,dt,z0,dir)
    r0=norm(R0); % length of R0
    r1=norm(R1); % length of R1

    dtheta=acos(dot(R0,R1)/r0/r1); %transfer angle

    if dir==1 %short way
        sgn=1; % "direction of motion"
    else
        dtheta=2*pi-dtheta; %long way
        sgn=-1;
    end

    A=sgn*sqrt(r0*r1*(1+cos(dtheta)));
    z=z0; %trial value
    relerr=1; %initialize error
    tol=1e-5; % Accepted relative error

    while relerr>tol
        if z==0 %to avoid division by zero
            S=1/6; Sprime=-1/120;
            C=1/2; Cprime=-1/24;
        else
            C=(1-cos(sqrt(z)))/z;
            S=(sqrt(z)-sin(sqrt(z)))/z^(3/2);
            Cprime=(z*sin(sqrt(z))-2*(sqrt(z)-1)*cos(sqrt(z)))/(2*z^(5/2));
            Sprime=(2*(1-cos(sqrt(z)))-z*sin(sqrt(z)))/(2*z^(5/2));
        end
        % calculate u(z), v(z), new value of z
        u=r0+r1+A*(z*S-1)/sqrt(C);
        uprime=A*(z*Sprime+Cprime)/sqrt(C);
        v=uprime;
        znew=z+(dt-sqrt(mu)*z^3*S-A*z^2*C)/v;
        relerr=abs(znew-z)/z;
        z=znew;


    end
    % use value of z obtained to determine v1,v2
    f=1-u/r0;
    gdot=1-u/r1;
    g=A*sqrt(u/mu);
    V1=(R1-f*R0)/g;
    V2=(gdot*R1-R0)/g;
    disp('V1=');
    disp(V1);
    disp('V2=');
    disp(V2);
end
 