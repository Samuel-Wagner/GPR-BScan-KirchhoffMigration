function [theta_a, theta_g, phi] = GPR_transmission_angles_v4(er,h,xa,ya,xq,yq,zq)
% find_GPR_transmission_angles
% is a speed > accuracy function to find the ray-based refraction
% angles for a GPR scenario. It is inaccurate with low depth.
% it uses an iterative Newton-Raphson method to solve for the angles
% which result in the light's "path of least time"


%er - relative permittivity of ground
%h  - height of antenna above ground
%d  - depth of buried target
%xA - antenna phase center x 
%yA - antenna phase center y
%xq - x position of target
%yq - y position of target
%zq - z position of target

%thetat - transmitted (in air) angle
%thetar - refracted (in soil) angle
%phi - azimuth angle

xa=single(xa);
ya=single(ya);
xq=single(xq);
yq=single(yq);

max_n_steps = 50;
yq=yq-ya;
xq=xq-xa;

phi=atan2(yq,xq);
if(phi < 0)
    phi = 2*pi + phi;
end
xi = xq; 
tol= 1e-4;
step=2*tol;

nsteps = 0;

while(abs(step)>tol) 

    %calculate value of F'(xi)
    t_center = xi*(tan(phi)^2 + 1)/sqrt(xi^2 + (xi*tan(phi))^2 + h^2) - ...
        sqrt(er)*(tan(phi)*(yq-xi*tan(phi)) + (xq-xi))/sqrt((xq-xi)^2+(yq-xi*tan(phi))^2 +zq^2);

    %calculate value of F'(xi+delta)
    t_right = (xi+tol)*(tan(phi)^2 + 1)/sqrt((xi+tol)^2 + ((xi+tol)*tan(phi))^2 + h^2) - ...
        sqrt(er)*(tan(phi)*(yq-(xi+tol)*tan(phi)) + (xq-(xi+tol)))/sqrt((xq-(xi+tol))^2+(yq-(xi+tol)*tan(phi))^2 +zq^2);

    %calculate step size as F'(xi)/F''(xi). /2 for stability
    step = -t_center/((t_right-t_center)/tol)/2;

    %update xi with step size
    xi = xi + step;

    nsteps = nsteps+1;
    if(nsteps > max_n_steps)
        xi = xq;
        break;
    end
end

theta_a=acos(h/sqrt(xi^2 + (xi.*tan(phi)).^2 + h^2));
theta_g=acos(zq/sqrt((xq-xi)^2 + (yq - xi.*tan(phi))^2 + zq^2));
