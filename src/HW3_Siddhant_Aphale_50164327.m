close all; clear; clc;
%Domain is a [-2, 2] x [-2, 2] with a 64 x 64 grid
xmax=2; xmin=-2; %Defining the maximum and minimum values of domain in X
ymax=2; ymin=-2; %Defining the maximum and minimum values of domain in y
imax=64; jmax=64; %Defining the grid
%Finding out the grid spacing values
dx=(xmax-xmin)/(imax-1);
dy=(ymax-ymin)/(jmax-1);
dt=0.5*dx; %Defining the time-step

%Finding out the values of X and Y at every grid points
[Y, X] = meshgrid(-2:dx:2, -2:dy:2);
% dlmwrite('X.txt',X);
% dlmwrite('Y.txt',Y);

%Total number of iterations
n=512;

%Initialization
for i=1:imax
    for j=1:jmax
%         phi_zero(i,j)=0.25.*(X(i,j).^2)+(Y(i,j).^2)-0.25;
        phi_zero(i,j)=16*(X(i,j).^4)+16*(Y(i,j).^4)-4*(X(i,j).^2)-4*(Y(i,j).^2);
%         phi_zero(i,j)=((X(i,j).^2)+(Y(i,j).^2)).^2-1;
    end
end
% dlmwrite('Initial.txt',phi_zero);
phi=phi_zero;
for t=1:n
    
    %Reinitialization of phi
    %Calculating the values of a,b,c and d first
    
    for i=2:imax-1
        for j=2:jmax-1
            a(i,j)=(phi(i,j)-phi(i-1,j))/dx;
            b(i,j)=(phi(i+1,j)-phi(i,j))/dx;
            c(i,j)=(phi(i,j)-phi(i,j-1))/dy;
            d(i,j)=(phi(i,j+1)-phi(i,j))/dy;
        end
    end
    
    %Calculation of a,b,c and d at boundaries by using Ghost Nodes
        
    for j=1:jmax
        a(1,j)=(phi(1,j)-phi(imax-1,j))/dx;
        b(1,j)=(phi(2,j)-phi(1,j))/dx;
        a(imax,j)=(phi(imax,j)-phi(imax-1,j))/dx;
        b(imax,j)=(phi(2,j)-phi(imax,j))/dx;
    end
    for i=1:imax
        c(i,1)=(phi(i,1)-phi(i,jmax-1))/dy;
        d(i,1)=(phi(i,2)-phi(i,1))/dy;
        c(i,jmax)=(phi(i,jmax)-phi(i,jmax-1))/dy;
        d(i,jmax)=(phi(i,2)-phi(i,jmax))/dy;
    end
    for i=2:imax-1
        a(i,1)=(phi(i,1)-phi(i-1,1))/dx;
        a(i,jmax)=(phi(i,jmax)-phi(i-1,jmax))/dx;
        b(i,1)=(phi(i+1,1)-phi(i,1))/dx;
        b(i,jmax)=(phi(i+1,jmax)-phi(i,jmax))/dx;
    end
    for j=2:jmax-1
        c(1,j)=(phi(1,j)-phi(1,j-1))/dy;
        c(imax,j)=(phi(imax,j)-phi(imax,j-1))/dy;
        d(1,j)=(phi(1,j+1)-phi(1,j))/dy;
        d(imax,j)=(phi(imax,j+1)-phi(imax,j))/dy;
    end
    
    %Finding out the values of a+, a-, b+, b-, c+, c-, d+, d-, s+, s-
    
    for i=1:imax
        for j=1:jmax
            a_plus(i,j)=max(a(i,j),0);
            a_minus(i,j)=min(0,a(i,j));
            b_plus(i,j)=max(b(i,j),0);
            b_minus(i,j)=min(0,b(i,j));
            c_plus(i,j)=max(c(i,j),0);
            c_minus(i,j)=min(0,c(i,j));
            d_plus(i,j)=max(d(i,j),0);
            d_minus(i,j)=min(0,d(i,j));
            
            if phi(i,j)>0.0001
                s_plus(i,j)=+1;
                s_minus(i,j)=0;
            elseif phi(i,j)<-0.0001
                s_plus(i,j)=0;
                s_minus(i,j)=-1;
            else
                s_plus(i,j)=phi(i,j)/(sqrt(phi(i,j)^2+dx^2));
                s_minus(i,j)=phi(i,j)/(sqrt(phi(i,j)^2+dx^2));
            end
            if phi(i,j)>0.0001
                s_zero(i,j)=1;
            end
            if phi(i,j)<-0.0001
                s_zero(i,j)=-1;
            end
            if phi(i,j)==0
                s_zero(i,j)=phi(i,j)/(sqrt(phi(i,j)^2+dx^2));
            end
        end
    end
     
    %Finding out the new values of phi
    
    for i=1:imax
        for j=1:jmax
            phi_new(i,j)=phi(i,j)-dt*(s_plus(i,j)*...
                ((sqrt(max(a_plus(i,j)^2,b_minus(i,j)^2)+...
                max(c_plus(i,j)^2,d_minus(i,j)^2)))-1))-...
                dt*(s_minus(i,j)*((sqrt(max(a_minus(i,j)^2,b_plus(i,j)^2)+...
                max(c_minus(i,j)^2,d_plus(i,j)^2)))-1));
        end
    end
   
    %Sub-cell Fix
    %Calculating signed distance functions
    for i=2:imax-1
        for j=2:jmax-1
            if phi(i+1,j)*phi(i,j)<0 ||...
                    phi(i,j+1)*phi(i,j)<0 ||...
                    phi(i,j)*phi(i-1,j)<0 ||...
                    phi(i,j)*phi(i,j-1)<0
                SDF(i,j)=(2*dx*phi_zero(i,j))/...
                    (sqrt((phi_zero(i+1,j)-phi_zero(i-1,j)).^2+...
                    (phi_zero(i,j+1)-phi_zero(i,j-1)).^2));
            else
                SDF(i,j)=0;
            end
        end
    end
    for i=2:imax-1
        if phi(i+1,1)*phi(i,1)<0||...
                phi(i,2)*phi(i,1)<0||...
                phi(i,1)*phi(i-1,1)<0||...
                phi(i,1)*phi(i,jmax-1)<0
            SDF(i,1)=(2*dx*phi_zero(i,1))/...
                (sqrt((phi_zero(i+1,1)-phi_zero(i-1,1)).^2+...
                (phi_zero(i,2)-phi_zero(i,jmax-1)).^2));
            phi_sub_fix(i,1)=phi(i,1)-dt*((s_zero(i,1)*...
                abs(phi(i,1)))-SDF(i,1))/dx;
        else
            SDF(i,1)=0;
            phi_sub_fix(i,1)=phi(i,1)-dt*(s_plus(i,1)*...
                ((sqrt(max(a_plus(i,1)^2,b_minus(i,1)^2)+...
                max(c_plus(i,1)^2,d_minus(i,1)^2)))-1))-...
                dt*(s_minus(i,1)*((sqrt(max(a_minus(i,1)^2,...
                b_plus(i,1)^2)+max(c_minus(i,1).^2,d_plus(i,1).^2)))-1));
        end
        if phi(i+1,jmax)*phi(i,jmax)<0||...
                phi(i,2)*phi(i,jmax)<0||...
                phi(i,jmax)*phi(i-1,jmax)<0||...
                phi(i,jmax)*phi(i,jmax-1)<0
            SDF(i,jmax)=(2*dx*phi_zero(i,jmax))/...
                (sqrt((phi_zero(i+1,jmax)-phi_zero(i-1,jmax)).^2+...
                (phi_zero(i,jmax-1)-phi_zero(i,2)).^2));
            phi_sub_fix(i,jmax)=phi(i,jmax)-dt.*((s_zero(i,jmax)*...
                abs(phi(i,jmax)))-SDF(i,jmax))/dx;
        else
            SDF(i,jmax)=0;
            phi_sub_fix(i,jmax)=phi(i,jmax)-dt*(s_plus(i,jmax)*...
                ((sqrt(max(a_plus(i,jmax)^2,b_minus(i,jmax)^2)+...
                max(c_plus(i,jmax)^2,d_minus(i,jmax)^2)))-1))-...
                dt*(s_minus(i,jmax)*((sqrt(max(a_minus(i,jmax)^2,...
                b_plus(i,jmax)^2)+max(c_minus(i,jmax)^2,d_plus(i,jmax)^2)))-1));
        end
    end
    
    for j=2:jmax-1
        if phi(2,j)*phi(1,j)<0||...
                phi(1,j+1)*phi(1,j)<0||...
                phi(1,j)*phi(imax-1,j)<0||...
                phi(1,j)*phi(1,j-1)<0
            SDF(1,j)=(2*dx*phi_zero(1,j))/...
                (sqrt((phi_zero(2,j)-phi_zero(imax-1,j)).^2+...
                ((phi-zero(1,j+1)-phi_zero(1,j-1)).^2)));
            phi_sub_fix(1,j)=phi(1,j)-dt*((s_zero(1,j)*...
                abs(phi(1,j)))-SDF(1,j))/dx;
        else
            SDF(1,j)=0;
            phi_sub_fix(1,j)=phi(1,j)-dt*(s_plus(1,j)*...
                ((sqrt(max(a_plus(1,j)^2,b_minus(1,j)^2)+...
                max(c_plus(1,j)^2,d_minus(1,j)^2)))-1))-...
                dt*(s_minus(1,j)*((sqrt(max(a_minus(1,j)^2,...
                b_plus(1,j)^2)+max(c_minus(1,j)^2,d_plus(1,j)^2)))-1));
        end
            if phi(2,j)*phi(imax,j)<0||...
                    phi(imax,j+1)*phi(imax,j)<0||...
                    phi(imax,j)*phi(imax-1,j)<0||...
                    phi(imax,j)*phi(imax,j-1)<0
            SDF(imax,j)=(2*dx*phi_zero(imax,j))/...
                (sqrt((phi_zero(2,j)-phi_zero(imax-1,j)).^2+...
                (phi_zero(imax,j+1)-phi_zero(imax,j-1)).^2));
            phi_sub_fix(imax,j)=phi(imax,j)-dt*((s_zero(imax,j)*...
                abs(phi(imax,j)))-SDF(imax,j))/dx;
        else
            SDF(imax,j)=0;
            phi_sub_fix(imax,j)=phi(imax,j)-dt*(s_plus(imax,j)*...
                ((sqrt(max(a_plus(imax,j)^2,b_minus(imax,j)^2)+...
                max(c_plus(imax,j)^2,d_minus(imax,j)^2)))-1))-...
                dt*(s_minus(imax,j)*((sqrt(max(a_minus(imax,j)^2,...
                b_plus(imax,j)^2)+max(c_minus(imax,j)^2,d_plus(imax,j)^2)))-1));
        end
    end
    
    %Calculating the value of phi after sub-cell fix by marching in time
    for i=2:imax-1
        for j=2:jmax-1
            if phi(i+1,j)*phi(i,j)<0||...
                    phi(i,j+1)*phi(i,j)<0||...
                    phi(i,j)*phi(i-1,j)<0||...
                    phi(i,j)*phi(i,j-1)<0
                phi_sub_fix(i,j)=phi(i,j)-dt*((s_zero(i,j)*...
                    abs(phi(i,j)))-SDF(i,j))/dx;
            else
                phi_sub_fix(i,j)=phi(i,j)-dt*(s_plus(i,j)*...
                    ((sqrt(max(a_plus(i,j)^2,b_minus(i,j)^2)+...
                    max(c_plus(i,j)^2,d_minus(i,j)^2)))-1))-dt*...
                    (s_minus(i,j)*((sqrt(max(a_minus(i,j)^2,...
                    b_plus(i,j)^2)+max(c_minus(i,j)^2,d_plus(i,j)^2)))-1));
            end
        end
    end
    phi_sub_fix(1,1)=phi_sub_fix(2,1);
    phi_sub_fix(imax,jmax)=phi_sub_fix(imax-1,jmax);
    phi_sub_fix(1,jmax)=phi_sub_fix(2,jmax);
    phi_sub_fix(imax,1)=phi_sub_fix(imax-1,1);
    
    %For sub-cell fix un-comment this
    phi=phi_sub_fix;
    %For without sub-cell fix un-comment this
%     phi=phi_new;
end
% dlmwrite('withoutsub.txt',phi);
%Plotting the Initial Condition
figure(1)
contour(X,Y,phi_zero)
grid on
axis square
hold on;
contour(X,Y,phi_zero,[0,0],'-k','LineWidth',2)

%Plotting the Reinitialization
figure(2)
contour(X,Y,phi,[0,0],'-b','LineWidth',2)
grid on
axis square
hold on;
% contour(X,Y,phi_zero,[0,0],'-r','LineWidth',2)
contour(X,Y,phi)