clc; clear; close all;
% Symplectic Scheme, 2D Boussinesq eqn, Composition method
% ASrinivasan, 13Oct2024, uses Laplacian operator

addpath('C:\Anand\Acer_Data\SDSU\MOLE\mole-master\mole_MATLAB')
Mult = 2; kX = 60;
m = Mult*64; xmin = -kX; xmax = kX; 
n = Mult*64; ymin = -kX; ymax = kX; 
dx = (xmax-xmin)/m; dy = (ymax-ymin)/n;
tEnd = 20; cV = 1; 
dt = 1; %CFL*dx/cV; 

% initial condition
xGrid = [xmin xmin+dx/2:dx:xmax-dx/2 xmax]';
yGrid = xGrid;
[X, Y] = meshgrid(xGrid, yGrid); 

CaseNo = 2; % only 2 works
[u0, v0] = u0Func(X, Y, CaseNo); 
u0 = [u0(:); v0(:)]; 

[uCOMP4, ECOMP4, tCOMP4] = COMP4(m,dx,n,dy,dt,tEnd,u0,cV); 
[uRRK, tRRK, gamRRK, eRRK] = RRK(0,m,dx,n,dy,dt,tEnd,u0,cV); 
dt = 0.2; 
[uCOMP6, ECOMP6, tCOMP6] = COMP61(m,dx,n,dy,dt,tEnd,u0,cV); 


% uOut1 = uTimeOut(uRRK, m, n, 1); uOut2 = uTimeOut(uRRK, m, n, 50); 
% uOut3 = uTimeOut(uRRK, m, n, 100); uOut4 = uTimeOut(uRRK, m, n, 250); 

uOut1 = uTimeOut(uCOMP4, m, n, 1); uOut2 = uTimeOut(uCOMP4, m, n, 10); 
% uOut3 = uTimeOut(uCOMP4, m, n, 50); uOut4 = uTimeOut(uCOMP4, m, n, 100); 

% uOut1 = uTimeOut(uCOMP6, m, n, 1); uOut2 = uTimeOut(uCOMP6, m, n, 50); 
% uOut3 = uTimeOut(uCOMP6, m, n, 100); uOut4 = uTimeOut(uCOMP6, m, n, 150); 


figure; 
tcl = tiledlayout(2,1); 
nexttile
    mesh(X, Y, uOut1, 'edgecolor', 'k');     
    title('t = 0 s'); 
    xlabel('X'); ylabel('Y'); zlabel('u')
nexttile
    mesh(X, Y, uOut2, 'edgecolor', 'k');     
    title('t = 10 s'); 
    xlabel('X'); ylabel('Y'); zlabel('u')
% nexttile
%     mesh(X, Y, uOut3, 'edgecolor', 'k');     
%     title('t = 20 s'); 
%     xlabel('X'); ylabel('Y'); zlabel('u')
% nexttile
%     mesh(X, Y, uOut4, 'edgecolor', 'k');     
%     title('t = 30 s'); 
%     xlabel('X'); ylabel('Y'); zlabel('u')
% title(tcl,'2D Wave Eq, MIM6-COMP6')


E0ref = eRRK(1); 
E0refCOMP = ECOMP4(1); 
E0refCOMP6 = ECOMP6(1); 
pn1 = 1; pn2 = 1; figure; 
plot(tCOMP4(pn1:pn2:end), abs(ECOMP4(pn1:pn2:end)/E0refCOMP), '-.k', 'LineWidth', 1.2);  %
hold on; plot(tCOMP6(pn1:3:end), abs(ECOMP6(pn1:3:end)/E0refCOMP6), 'sr', 'LineWidth', 1.2); 
hold on; plot(tRRK(pn1:pn2:end), abs(eRRK(pn1:pn2:end)/E0ref), '--b', 'LineWidth', 1.2); %
title('Energy En/E0 vs time, 2D Boussinesq Eq'); % set(gca, 'yscale', 'log');
xlabel('time(s)'); ylabel('En/E0'); 
legend('MIM4-COMP4', 'MIM6-COMP6', 'MIM4-RK4', 'location', 'best')
set(gca,'FontSize',10); grid on; %ylim([0.9960, 1.0005]);



function [uCOMP4, ECOMP4, tCOMP4] = COMP4(m,dx,n,dy,dt,tEnd,u0, cV)

% 4th order Composition method
% Hairer et al, Geometric Integration, pp152, Sec V.3, eq. (3.6)

fprintf('Executing COMP4 \n'); 
U = u0(1:(m+2)*(n+2),1); Ut = u0((m+2)*(n+2)+1:end,1); 
V = 0*U; 
uCOMP4(:,1) = [U; Ut]; % initial value

    D = div2D(4, m, dx, n, dy);
    G = grad2D(4, m, dx, n, dy);
    L = cV^2*D*G; 
    ID = interpDMat2D(4,m,n); IG = interpGMat2D(4,m,n);
    D = D*ID; G = IG*G; 
    I1 = speye((m+2)*(n+2), (m+2)*(n+2));


alp1 = (146 + 5*sqrt(19))/540;  bet5 = alp1;
alp2 = (-2 + 10*sqrt(19))/135;  bet4 = alp2; 
alp3 = 1/5;                     bet3 = alp3;
alp4 = (-23 - 20*sqrt(19))/270; bet2 = alp4;
alp5 = (14 - sqrt(19))/108;     bet1 = alp5; 

qx = V; qxx = G*V; 
E1 = 0.5*(qx'*qx); E2 = 0.5*(qxx'*qxx); 
E3 = 0.5*(U'*U); E4 = 1/3*(U.^2)'*U; 

E0 = dx*dy*(E1 + E2 + E3 + E4);

ECOMP4(1,1) = E0; 
tCOMP4(1,1) = 0; % initial time    
iCount = 2; 
t = dt; 
V = Ut; 

    while t <= tEnd
        p1 = U + bet1*dt*G*V; 
        q1 = V + (I1 - L)\((bet1 + alp1)*dt*(G*(p1 + p1.^2)));         
        
        % Step 2        
        p2 = p1 + (bet2 + alp1)*dt*G*q1; 
        q2 = q1 + (I1 - L)\((bet2 + alp2)*dt*(G*(p2 + p2.^2))); 
    
    %
        % Step 3    
        p3 = p2 + (bet3 + alp2)*dt*G*q2; 
        q3 = q2 + (I1 - L)\((bet3 + alp3)*dt*(G*(p3 + p3.^2))); 
    
    
        % Step 4    
        p4 = p3 + (bet4 + alp3)*dt*G*q3; 
        q4 = q3 + (I1 - L)\((bet4 + alp4)*dt*(G*(p4 + p4.^2))); 

        % Step 5    
        p5 = p4 + (bet5 + alp4)*dt*G*q4; 
        q5 = q4 + (I1 - L)\((bet5 + alp5)*dt*(G*(p5 + p5.^2))); 
        
        p5 = p5 + alp5*dt*G*q5; 

        uNew = [p5; q5];                          
        
        qx = q5; qxx = G*q5; 
        E1 = 0.5*(qx'*qx); E2 = 0.5*(qxx'*qxx); 
        E3 = 0.5*(p5'*p5); E4 = 1/3*(p5.^2)'*p5; 

        ECOMP4(iCount,1) = dx*dy*(E1 + E2 + E3 + E4) ; 

        uCOMP4(:,iCount) = uNew;   
        U = uCOMP4(1:(m+2)*(n+2), iCount);
        V = uCOMP4((m+2)*(n+2)+1:end, iCount);

        tCOMP4(iCount,1) = t;
        
        iCount = iCount+1; 
        t = t + dt;              
        
    end

end


function [uCOMP4, ECOMP4, tCOMP4] = COMP61(m,dx,n,dy,dt,tEnd,u0, cV)

% 6th order Composition Method
% Hairer et al, Geometric Integration, pp152, Sec V.3
% Eq (3.9), (3.11)

fprintf('Executing COMP6 \n'); 
U = u0(1:(m+2)*(n+2),1); Ut = u0((m+2)*(n+2)+1:end,1); 
V = 0*U; 
uCOMP4(:,1) = [U; Ut]; % initial value

    D = div2D(6, m, dx, n, dy);
    G = grad2D(6, m, dx, n, dy);
    L = cV^2*D*G; 
    ID = interpDMat2D(6,m,n); IG = interpGMat2D(6,m,n);
    D = D*ID; G = IG*G; 
    I1 = speye((m+2)*(n+2), (m+2)*(n+2));

gam1 = 0.78451361047755726381949763;
gam2 = 0.23557321335935813368479318;
gam3 = -1.17767998417887100694641568;
gam4 = 1.31518632068391121888424973;
gam5 = gam3; gam6 = gam2; gam7 = gam1; 


alp1 = gam1/2;  bet1 = gam1/2; 
alp2 = gam2/2;  bet2 = gam2/2; 
alp3 = gam3/2;  bet3 = gam3/2; 
alp4 = gam4/2;  bet4 = gam4/2; 
alp5 = gam5/2;  bet5 = gam5/2; 
alp6 = gam6/2;  bet6 = gam6/2; 
alp7 = gam7/2;  bet7 = gam7/2; 


qx = V; qxx = G*V; 
E1 = 0.5*(qx'*qx); E2 = 0.5*(qxx'*qxx); 
E3 = 0.5*(U'*U); E4 = 1/3*(U.^2)'*U; 

E0 = dx*dy*(E1 + E2 + E3 + E4);
ECOMP4(1,1) = E0; 
tCOMP4(1,1) = 0; % initial time    
iCount = 2; 
t = dt; V = Ut; 

    while t <= tEnd
        p1 = U + bet1*dt*G*V; 
        q1 = V + (I1 - L)\((bet1 + alp1)*dt*(G*(p1 + p1.^2)));         
        
        % Step 2        
        p2 = p1 + (bet2 + alp1)*dt*G*q1; 
        q2 = q1 + (I1 - L)\((bet2 + alp2)*dt*(G*(p2 + p2.^2))); 
    
    %
        % Step 3    
        p3 = p2 + (bet3 + alp2)*dt*G*q2; 
        q3 = q2 + (I1 - L)\((bet3 + alp3)*dt*(G*(p3 + p3.^2))); 
    
    
        % Step 4    
        p4 = p3 + (bet4 + alp3)*dt*G*q3; 
        q4 = q3 + (I1 - L)\((bet4 + alp4)*dt*(G*(p4 + p4.^2))); 

        % Step 5    
        p5 = p4 + (bet5 + alp4)*dt*G*q4; 
        q5 = q4 + (I1 - L)\((bet5 + alp5)*dt*(G*(p5 + p5.^2))); 

        % Step 6    
        p6 = p5 + (bet6 + alp5)*dt*G*q5; 
        q6 = q5 + (I1 - L)\((bet6 + alp6)*dt*(G*(p6 + p6.^2)));

        % Step 7
        p7 = p6 + (bet7 + alp6)*dt*G*q6; 
        q7 = q6 + (I1 - L)\((bet7 + alp7)*dt*(G*(p7 + p7.^2)));

        p7 = p7 + alp7*dt*G*q7; 

        uNew = [p7; q7];                  
        
        qx = q7; qxx = G*q7; 
        E1 = 0.5*(qx'*qx); E2 = 0.5*(qxx'*qxx); 
        E3 = 0.5*(p5'*p5); E4 = 1/3*(p5.^2)'*p5; 

        ECOMP4(iCount,1) = dx*dy*(E1 + E2 + E3 + E4) ; 
        uCOMP4(:,iCount) = uNew;   
        U = uCOMP4(1:(m+2)*(n+2), iCount);
        V = uCOMP4((m+2)*(n+2)+1:end, iCount);

        tCOMP4(iCount,1) = t;
        
        iCount = iCount+1; 
        t = t + dt;              
        
    end

end



function [uCOMP4, ECOMP4, tCOMP4] = COMP62(m,dx,n,dy,dt,tEnd,u0, cV)

% 6th order Composition Method
% Hairer et al, Geometric Integration, pp152, Sec V.3
% Eq (3.9), (3.12)

U = u0(1:(m+2)*(n+2),1); V = u0((m+2)*(n+2)+1:end,1); 
uCOMP4(:,1) = [U; V]; % initial value

    D = div2D(6, m, dx, n, dy);
    G = grad2D(6, m, dx, n, dy);
    L = cV^2*D*G; 

gam1 = 0.39216144400731413927925056;
gam2 = 0.33259913678935943859974864;
gam3 = -0.70624617255763935980996482;
gam4 = 0.08221359629355080023149045;
gam5 = 0.79854399093482996339895035; 
gam6 = gam4; gam7 = gam3; 
gam8 = gam2; gam9 = gam1; 


alp1 = gam1/2;  bet1 = gam1/2; 
alp2 = gam2/2;  bet2 = gam2/2; 
alp3 = gam3/2;  bet3 = gam3/2; 
alp4 = gam4/2;  bet4 = gam4/2; 
alp5 = gam5/2;  bet5 = gam5/2; 
alp6 = gam6/2;  bet6 = gam6/2; 
alp7 = gam7/2;  bet7 = gam7/2; 
alp8 = gam8/2;  bet8 = gam8/2; 
alp9 = gam9/2;  bet9 = gam9/2; 


ut = cV*G*U;
E0 = 0.5*dx*dy*(ut'*ut + V'*V + 2*(1. - cos(U))'*ones(size(U,1), 1));
ECOMP4(1,1) = 0; 
tCOMP4(1,1) = 0; % initial time    
iCount = 2; 
t = dt; 

    while t <= tEnd
        p1 = U + bet1*dt*V; 
        q1 = V + (bet1 + alp1)*dt*(L*p1 - sin(p1));         
        
        % Step 2        
        p2 = p1 + (bet2 + alp1)*dt*q1; 
        q2 = q1 + (bet2 + alp2)*dt*(L*p2 - sin(p2)); 
    
    %
        % Step 3    
        p3 = p2 + (bet3 + alp2)*dt*q2; 
        q3 = q2 + (bet3 + alp3)*dt*(L*p3 - sin(p3)); 
    
    
        % Step 4    
        p4 = p3 + (bet4 + alp3)*dt*q3; 
        q4 = q3 + (bet4 + alp4)*dt*(L*p4 - sin(p4)); 

        % Step 5    
        p5 = p4 + (bet5 + alp4)*dt*q4; 
        q5 = q4 + (bet5 + alp5)*dt*(L*p5 - sin(p5));

        % Step 6    
        p6 = p5 + (bet6 + alp5)*dt*q5; 
        q6 = q5 + (bet6 + alp6)*dt*(L*p6 - sin(p6));

        % Step 7
        p7 = p6 + (bet7 + alp6)*dt*q6; 
        q7 = q6 + (bet7 + alp7)*dt*(L*p7 - sin(p7));

        % Step 8
        p8 = p7 + (bet8 + alp7)*dt*q7; 
        q8 = q7 + (bet8 + alp8)*dt*(L*p8 - sin(p8));

        % Step 9
        p9 = p8 + (bet9 + alp8)*dt*q8; 
        q9 = q8 + (bet9 + alp9)*dt*(L*p9 - sin(p9));


        p9 = p9 + alp9*dt*q9; 

        uNew = [p9; q9];                  
        
        ut1 = cV*G*p9; 
        ECOMP4(iCount,1) = 0.5*dx*dy*(ut1'*ut1 + q9'*q9 + 2*(1.- cos(p9))'*ones(size(p9,1), 1)) - E0; 
        uCOMP4(:,iCount) = uNew;   
        U = uCOMP4(1:(m+2)*(n+2), iCount);
        V = uCOMP4((m+2)*(n+2)+1:end, iCount);

        tCOMP4(iCount,1) = t;
        
        iCount = iCount+1; 
        t = t + dt;              
        
    end

end



function [uRRK, tRRK, gamRRK, eRRK] = RRK(RKFlag,m,dx,n,dy,dt,tEnd,Uin0, cV)
    
% Relaxation RK4
    fprintf('Executing RK4 \n'); 
    
    mn = (m+2)*(n+2); 
    D = div2D(4, m, dx, n, dy);
    G = grad2D(4, m, dx, n, dy);
    L = cV^2*D*G; 
    ID = interpDMat2D(4,m,n); IG = interpGMat2D(4,m,n);
    D = D*ID; G = IG*G; 
    I1 = speye(mn, mn);
              
    C(1,1) = 0;    C(2,1) = 1/2;
    C(3,1) = 1/2;    C(4,1) = 1;

    A(2,1) = 1/2;
    A(3,1) = 0;   A(3,2) = 1/2;
    A(4,1) = 0;   A(4,2) = 0;     A(4,3) = 1;     

    B(1,1) = 1/6;    B(2,1) = 1/3;
    B(3,1) = 1/3;    B(4,1) = 1/6;
     
    tRRK(1,1) = 0; % initial time    
    gamRRK(1,1) = 1;
    
    u0 = Uin0(1:mn,1); ut0 = Uin0(mn+1:end,1);     
    v0 = zeros(size(u0)); 

    qx = v0; qxx = G*v0;     
    E1 = 0.5*(qx'*qx); E2 = 0.5*(qxx'*qxx); 
    E3 = 0.5*(u0'*u0); E4 = 1/3*(u0.^2)'*u0; 
    E0 = dx*dy*(E1 + E2 + E3 + E4);
    
    eRRK(1,1) = E0; 
    uRRK(:,1) = [u0; ut0]; % initial value
    y = [u0; ut0];
    iCount = 2; 
    t = dt; 
   
    while t <= tEnd

        z1 = y; 
        z1U = z1(1:mn); z1V = z1(mn+1:end); 
        k1U = G*z1V; k1V = (I1 - L)\(G*(z1U + z1U.^2)); 
        k1 = [k1U; k1V]; %
        
        % Step 2        
        z2 =  y + dt*A(2,1)*k1;       
        z2U = z2(1:mn); z2V = z2(mn+1:end); 
        k2U = G*z2V; k2V = (I1 - L)\(G*(z2U + z2U.^2)); 
        k2 = [k2U; k2V]; %

        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
        z3U = z3(1:mn); z3V = z3(mn+1:end); 
        k3U = G*z3V; k3V = (I1 - L)\(G*(z3U + z3U.^2)); 
        k3 = [k3U; k3V]; %

        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        z4U = z4(1:mn); z4V = z4(mn+1:end); 
        k4U = G*z4V; k4V = (I1 - L)\(G*(z4U + z4U.^2)); 
        k4 = [k4U; k4V]; %
       
        BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
            B(4,1)*k4;    
    
        gam = 1;
        uNew = y + gam*dt*BKsum;  

        uRRK(:,iCount) = uNew;   
        y = uRRK(:,iCount);
        tRRK(iCount,1) = t;
        gamRRK(iCount,1) = gam; 
        
        Uin = uNew(1:mn,1); Vin = uNew(mn+1:end,1);
        qx = Vin; qxx = G*Vin; 
        E1 = 0.5*(qx'*qx); E2 = 0.5*(qxx'*qxx); 
        E3 = 0.5*(Uin'*Uin); E4 = 1/3*(Uin.^2)'*Uin; 

        eRRK(iCount,1) = dx*dy*(E1 + E2 + E3 + E4) ; 
        
        iCount = iCount+1; 
        t = t + gam*dt;           
        
    end   

end

    
function [u0, v0] = u0Func(X, Y, CaseNo)

switch CaseNo
    case 1
    % case 1 - refer Rigge , Minzoni
        lmb = 0.5; phi = 0.5; xi = 0; kk = 0;
        E1 = lmb/sqrt(1-lmb^2)*sin(phi - kk*X);
        E2 = sech(lmb*(X - xi));
        E3 = E2.*sech(lmb*Y); 
        
        u0 = -4*atan(E1.*E3);        
        v0 = zeros(size(u0)); 

    case 2
        
        A = 1; sigma = 2; mu = 0;         
        u0 = A*exp(-1/(sigma^2)*((Y-mu/2 ).^2 + (X-mu/2).^2));   
%         v0 = -u0; 
        v0 = zeros(size(u0)); 

    
end

end


function [yR] = tInterp(tR, uR, tEnd, NElem)

    t1 = tR(end-1, 1); t2 = tR(end, 1);
    u1 = uR(1:NElem+2, end-1); u2 = uR(1:NElem+2, end); 
    yR = u1 + (u2 - u1)*(tEnd - t1)/(t2 - t1); 

end

function uOut = uExact(t, xGrid, CaseNo)

switch CaseNo
    case 1
        uOut = 4*atan(t*sech(xGrid)); 
        
    case 2
        c = 0.5; B = (xGrid - c*t)/sqrt(1-c^2); 
        uOut = 4*atan(exp(B)); 

end

          

end

function uOut = uTimeOut(uRRK, m, n, tOut)

    uOut = uRRK(1:(m+2)*(n+2), tOut); 
    uOut = reshape(uOut, [m+2, n+2]); 


end