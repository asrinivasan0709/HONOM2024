clc; clear; close all;
% Symplectic Scheme, 3D Wave eqn, Composition method
% ASrinivasan, 6Apr2024

addpath('...\mole_MATLAB')
Mult = 1; 
m = Mult*64; n = m; o = m; 
xmin = -25; xmax = 25; 
ymin = -25; ymax = 25; zmin = -25; zmax = 25; 
dx = (xmax-xmin)/m; dy = (ymax-ymin)/n;  dz = (zmax-zmin)/o;
tEnd = 5; cV = 1; 
CFL = 0.10; dt = 0.1; % CFL*dx/cV; 

% initial condition
xGrid = [xmin xmin+dx/2:dx:xmax-dx/2 xmax]';
yGrid = xGrid; zGrid = xGrid; 
[X, Y, Z] = meshgrid(xGrid, yGrid, xGrid); 

CaseNo = 1; 
[u0, v0] = u0Func(X, Y, Z, CaseNo); 
u0 = [u0(:); v0(:)]; 

[uCOMP4, ECOMP4, tCOMP4] = COMP4(m,dx,n,dy,o,dz,dt,tEnd,u0,cV); 
[uCOMP6, ECOMP6, tCOMP6] = COMP61(m,dx,n,dy,o,dz,dt,tEnd,u0,cV); 
[uRRK, tRRK, gamRRK, eRRK] = RRK(0,m,dx,n,dy,o,dz,dt,tEnd,u0,cV); 

mno = (m+2)*(n+2)*(o+2); 
uEx = uExact(tRRK(end), X, Y, Z, CaseNo); uEx = uEx(:); 
uRRKOut = uRRK(1:mno, end); 
uCOMP4Out = uCOMP4(1:mno, end); 
uCOMP6Out = uCOMP6(1:mno, end); 

normOut = [m, dx, dt, norm(uEx - uRRKOut, 'inf'), ...
                norm(uEx - uCOMP4Out, 'inf'), ...
                norm(uEx - uCOMP6Out, 'inf')]; 

pn1 = 1; pn2 = 3; figure; 
plot(tCOMP4(pn1:pn2:end), abs(ECOMP4(pn1:pn2:end)/ECOMP4(1)), '-.k', 'LineWidth', 1.2); 
hold on; plot(tCOMP6(pn1:pn2:end), abs(ECOMP6(pn1:pn2:end)/ECOMP6(1)), 'sr', 'LineWidth', 1.2); 
hold on; plot(tRRK(pn1:pn2:end), abs(eRRK(pn1:pn2:end)/eRRK(1)), '--b', 'LineWidth', 1.2); 
title('Energy En/E0 vs time, Wave Eq 3D'); % set(gca, 'yscale', 'log');
xlabel('time(s)'); ylabel('En/E0'); legend('MIM4-COMP4', 'MIM6-COMP6', 'MIM4-RK4', 'location', 'best')
set(gca,'FontSize',10); grid on; %ylim([0.99995, 1.00001]); 


function [uCOMP4, ECOMP4, tCOMP4] = COMP4(m,dx,n,dy,o,dz,dt,tEnd,u0, cV)

% 4th order Composition method
% Hairer et al, Geometric Integration, pp152, Sec V.3, eq. (3.6)
fprintf('Executing COMP4 ... \n'); 
mno = (m+2)*(n+2)*(o+2); 

U = u0(1:mno,1); V = u0(mno+1:end,1); 
uCOMP4(:,1) = [U; V]; % initial value

    D = div3D(4, m, dx, n, dy, o, dz);
    G = grad3D(4, m, dx, n, dy, o, dz);
    ID = interpDMat3D(4,m,n,o); 
    IG = interpGMat3D(4,m,n,o); 
    D = D*ID; G = IG*G; 


alp1 = (146 + 5*sqrt(19))/540;  bet5 = alp1;
alp2 = (-2 + 10*sqrt(19))/135;  bet4 = alp2; 
alp3 = 1/5;                     bet3 = alp3;
alp4 = (-23 - 20*sqrt(19))/270; bet2 = alp4;
alp5 = (14 - sqrt(19))/108;     bet1 = alp5; 

E0 = 0.5*dx*dy*dz*(U'*U + V'*V );
ECOMP4(1,1) = E0; 
tCOMP4(1,1) = 0; % initial time    
iCount = 2; 
t = dt; 


    while t <= tEnd
        p1 = U + bet1*dt*D*V; 
        q1 = V + (bet1 + alp1)*dt*(G*p1 );         
        
        % Step 2        
        p2 = p1 + (bet2 + alp1)*dt*D*q1; 
        q2 = q1 + (bet2 + alp2)*dt*(G*p2 ); 
    
    %
        % Step 3    
        p3 = p2 + (bet3 + alp2)*dt*D*q2; 
        q3 = q2 + (bet3 + alp3)*dt*(G*p3 ); 
    
    
        % Step 4    
        p4 = p3 + (bet4 + alp3)*dt*D*q3; 
        q4 = q3 + (bet4 + alp4)*dt*(G*p4 ); 

        % Step 5    
        p5 = p4 + (bet5 + alp4)*dt*D*q4; 
        q5 = q4 + (bet5 + alp5)*dt*(G*p5 );
        
        p5 = p5 + alp5*dt*D*q5; 

        uNew = [p5; q5];                  
        
        ECOMP4(iCount,1) = 0.5*dx*dy*dz*(p5'*p5 + q5'*q5 ); % - E0; 
        uCOMP4(:,iCount) = uNew;   
        U = uCOMP4(1:mno, iCount);
        V = uCOMP4(mno+1:end, iCount);

        tCOMP4(iCount,1) = t;
        
        iCount = iCount+1; 
        t = t + dt;              
        
    end

end


function [uCOMP4, ECOMP4, tCOMP4] = COMP61(m,dx,n,dy,o,dz,dt,tEnd,u0, cV)

% 6th order Composition Method
% Hairer et al, Geometric Integration, pp152, Sec V.3
% Eq (3.9), (3.11)
fprintf('Executing COMP6 ...\n'); 
mno = (m+2)*(n+2)*(o+2);
U = u0(1:mno,1); V = u0(mno+1:end,1); 
uCOMP4(:,1) = [U; V]; % initial value

    D = div3D(6, m, dx, n, dy, o, dz);
    G = grad3D(6, m, dx, n, dy, o, dz);
    ID = interpDMat3D(6,m,n,o); 
    IG = interpGMat3D(6,m,n,o); 
    D = D*ID; G = IG*G; 

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


E0 = 0.5*dx*dy*dz*(U'*U + V'*V );
ECOMP4(1,1) = E0; 
tCOMP4(1,1) = 0; % initial time    
iCount = 2; 
t = dt; 

    while t <= tEnd
        p1 = U + bet1*dt*D*V; 
        q1 = V + (bet1 + alp1)*dt*(G*p1 );         
        
        % Step 2        
        p2 = p1 + (bet2 + alp1)*dt*D*q1; 
        q2 = q1 + (bet2 + alp2)*dt*(G*p2 ); 
    
    %
        % Step 3    
        p3 = p2 + (bet3 + alp2)*dt*D*q2; 
        q3 = q2 + (bet3 + alp3)*dt*(G*p3 ); 
    
    
        % Step 4    
        p4 = p3 + (bet4 + alp3)*dt*D*q3; 
        q4 = q3 + (bet4 + alp4)*dt*(G*p4 ); 

        % Step 5    
        p5 = p4 + (bet5 + alp4)*dt*D*q4; 
        q5 = q4 + (bet5 + alp5)*dt*(G*p5 );

        % Step 6    
        p6 = p5 + (bet6 + alp5)*dt*D*q5; 
        q6 = q5 + (bet6 + alp6)*dt*(G*p6 );

        % Step 7
        p7 = p6 + (bet7 + alp6)*dt*D*q6; 
        q7 = q6 + (bet7 + alp7)*dt*(G*p7 );

        p7 = p7 + alp7*dt*D*q7; 

        uNew = [p7; q7];                  
        

        ECOMP4(iCount,1) = 0.5*dx*dy*dz*(p7'*p7 + q7'*q7 ); % - E0; 
        uCOMP4(:,iCount) = uNew;   
        U = uCOMP4(1:mno, iCount);
        V = uCOMP4(mno+1:end, iCount);

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
E0 = 0.5*dx*dy*(ut'*ut + V'*V );
ECOMP4(1,1) = 0; 
tCOMP4(1,1) = 0; % initial time    
iCount = 2; 
t = dt; 

    while t <= tEnd
        p1 = U + bet1*dt*V; 
        q1 = V + (bet1 + alp1)*dt*(L*p1);         
        
        % Step 2        
        p2 = p1 + (bet2 + alp1)*dt*q1; 
        q2 = q1 + (bet2 + alp2)*dt*(L*p2 ); 
    
    %
        % Step 3    
        p3 = p2 + (bet3 + alp2)*dt*q2; 
        q3 = q2 + (bet3 + alp3)*dt*(L*p3 ); 
    
    
        % Step 4    
        p4 = p3 + (bet4 + alp3)*dt*q3; 
        q4 = q3 + (bet4 + alp4)*dt*(L*p4 ); 

        % Step 5    
        p5 = p4 + (bet5 + alp4)*dt*q4; 
        q5 = q4 + (bet5 + alp5)*dt*(L*p5 );

        % Step 6    
        p6 = p5 + (bet6 + alp5)*dt*q5; 
        q6 = q5 + (bet6 + alp6)*dt*(L*p6 );

        % Step 7
        p7 = p6 + (bet7 + alp6)*dt*q6; 
        q7 = q6 + (bet7 + alp7)*dt*(L*p7 );

        % Step 8
        p8 = p7 + (bet8 + alp7)*dt*q7; 
        q8 = q7 + (bet8 + alp8)*dt*(L*p8 );

        % Step 9
        p9 = p8 + (bet9 + alp8)*dt*q8; 
        q9 = q8 + (bet9 + alp9)*dt*(L*p9 );


        p9 = p9 + alp9*dt*q9; 

        uNew = [p9; q9];                  
        
        ut1 = cV*G*p9; 
        ECOMP4(iCount,1) = 0.5*dx*dy*(ut1'*ut1 + q9'*q9 ) - E0; 
        uCOMP4(:,iCount) = uNew;   
        U = uCOMP4(1:(m+2)*(n+2), iCount);
        V = uCOMP4((m+2)*(n+2)+1:end, iCount);

        tCOMP4(iCount,1) = t;
        
        iCount = iCount+1; 
        t = t + dt;              
        
    end

end



function [uRRK, tRRK, gamRRK, eRRK] = RRK(RKFlag,m,dx,n,dy,o,dz,dt,tEnd,u0, cV)
    
% Relaxation RK4
fprintf('Executing RK4 ... \n'); 

    mno = (m+2)*(n+2)*(o+2);
     
    D = div3D(4, m, dx, n, dy, o, dz);
    G = grad3D(4, m, dx, n, dy, o, dz);
    ID = interpDMat3D(4,m,n,o); 
    IG = interpGMat3D(4,m,n,o); 
    D = D*ID; G = IG*G; 
              
    C(1,1) = 0;    C(2,1) = 1/2;
    C(3,1) = 1/2;    C(4,1) = 1;

    A(2,1) = 1/2;
    A(3,1) = 0;   A(3,2) = 1/2;
    A(4,1) = 0;   A(4,2) = 0;     A(4,3) = 1;     

    B(1,1) = 1/6;    B(2,1) = 1/3;
    B(3,1) = 1/3;    B(4,1) = 1/6;
     
    uRRK(:,1) = u0; % initial value
    y = u0; 
    tRRK(1,1) = 0; % initial time    
    gamRRK(1,1) = 1;
    
    Uin0 = y(1:mno,1); Vin0 = y(mno+1:end,1); 
    E0 = 0.5*dx*dy*dz*(Uin0'*Uin0 + ...
            Vin0'*Vin0 ); 
    eRRK(1,1) = E0; 
    
    iCount = 2; 
    t = dt; 
   
    while t <= tEnd

        z1 = y; 
        z1U = z1(1:mno); z1V = z1(mno+1:end); 
        [k1] = [D*z1V; 
                G*z1U ]; %
        
        % Step 2        
        z2 =  y + dt*A(2,1)*k1;       
        z2U = z2(1:mno); z2V = z2(mno+1:end); 
        [k2] = [D*z2V; 
                G*z2U ]; %

        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
        z3U = z3(1:mno); z3V = z3(mno+1:end); 
        [k3] = [D*z3V; 
                G*z3U ]; %      

        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        z4U = z4(1:mno); z4V = z4(mno+1:end); 
        [k4] = [D*z4V; 
                G*z4U ]; %   

       
        BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
            B(4,1)*k4; % + B(5,1)*k5;   
    
        switch RKFlag
            case 1
                U = cV*y(1:NElem+2,1); V = y(NElem+3:end,1);
                dU = cV*BKsum(1:NElem+2,1); dV = BKsum(NElem+3:end,1);
                
                AA = U'*G'*G*U + V'*V;
                BB = 2*(U'*G'*G*dU + V'*dV);
                CC = dU'*G'*G*dU + dV'*dV;
                
                gam = (E0 - AA - BB)/(dt* CC);

            case 0
                gam = 1;
        end
        uNew = y + gam*dt*BKsum;  

        uRRK(:,iCount) = uNew;   
        y = uRRK(:,iCount);
        tRRK(iCount,1) = t;
        gamRRK(iCount,1) = gam; 
        
        Uin = uNew(1:mno,1); Vin = uNew(mno+1:end,1);
       eRRK(iCount,1) = 0.5*dx*dy*dz*(Uin'*Uin + ...
            Vin'*Vin ); % - E0;  
        
        iCount = iCount+1; 
        t = t + gam*dt;           
        
    end   

end






function [u0, v0] = u0Func(X, Y, Z, CaseNo)

switch CaseNo
    case 1
    % case 1 - refer Rigge 
        A = 1; sigma = 0.5; mu = 0;         
        u0 = A*exp(-1/(sigma^2)*((Z-mu/2 ).^2 + ...
                (Y-mu/2 ).^2 + (X-mu/2).^2));   
        v0 = -u0; 

    
end

end


function [yR] = tInterp(tR, uR, tEnd, NElem)

    t1 = tR(end-1, 1); t2 = tR(end, 1);
    u1 = uR(1:NElem+2, end-1); u2 = uR(1:NElem+2, end); 
    yR = u1 + (u2 - u1)*(tEnd - t1)/(t2 - t1); 

end

function uOut = uExact(t, X, Y, Z, CaseNo)

switch CaseNo
    case 1
        A = 1; sigma = 0.5; mu = 0;         
        uOut = A*exp(-1/(sigma^2)*((Z-mu/2 - t).^2 + ... 
                (Y-mu/2 - t).^2 + (X-mu/2 - t).^2)); 

end

          

end

function uOut = uTimeOut(uRRK, m, n, tOut)

    uOut = uRRK(1:(m+2)*(n+2), tOut); 
    uOut = reshape(uOut, [m+2, n+2]); 


end
