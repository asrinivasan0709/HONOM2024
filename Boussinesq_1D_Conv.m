clc; clear; close all;
% Symplectic Scheme, 1D Improved Boussinesq eqn, Composition method
% ASrinivasan, 4Apr2024
% Ref: Li et al, J. Comp Physics 2020
    % NOTE - use this file for convergence calcs

addpath('...\mole-master\mole_MATLAB')

mVec = [1,2,4,8]';
for i = 1:size(mVec,1)
    fprintf('Executing iteration %d of %d ... \n', i, size(mVec,1)); 

    Mult = mVec(i); 
    m = Mult*64; xmin = -50; xmax = 50; 
    dx = (xmax-xmin)/m; tEnd = 5; cV = 1; 
    CFL = 0.5; dt =  CFL*dx/cV; 
    
    % initial condition
    xGrid = [xmin xmin+dx/2:dx:xmax-dx/2 xmax]';
    xNod = [xmin:dx:xmax]';
    CaseNo = 1; % 2 works
    [u0, ut0] = u0Func(xGrid, CaseNo);  
    u0 = [u0; ut0]; 
    
    [uCOMP4, ECOMP4, tCOMP4] = COMP4(m,dx,dt,tEnd,u0,cV); 
    [uRRK, tRRK, gamRRK, eRRK] = RRK(0,m,dx,dt,tEnd,u0,cV); 
    % dt = 0.25;
    [uCOMP6, ECOMP6, tCOMP6] = COMP61(m,dx,dt,tEnd,u0,cV); 
    [uEx, utEx] = uExact(tRRK(end), xGrid, CaseNo); 
    
    normCalc = [norm(uRRK(1:m+2, end) - uEx, 'inf'), ...
                    norm(uCOMP4(1:m+2, end) - uEx, 'inf'), ...
                    norm(uCOMP6(1:m+2, end) - uEx, 'inf')]; 
    
    normOut(i, :) = [m, dx, dt, normCalc]; 
end

% Calculate order of convergence
for i = 1:size(mVec, 1) - 1
    CvRK4(i, 1) = 1/log(2)*log(normOut(i,4)/normOut(i+1,4)); 
    CvCOMP4(i, 1) = 1/log(2)*log(normOut(i,5)/normOut(i+1,5)); 
    CvCOMP6(i, 1) = 1/log(2)*log(normOut(i,6)/normOut(i+1,6)); 
    
end

Cv = [[0; CvRK4], [0; CvCOMP4], [0; CvCOMP6]];
normOut = [normOut, Cv]; 


figure;
loglog(normOut(:, 2), normOut(:, 5), '-*k');
hold on; loglog(normOut(:, 2), normOut(:, 6), '-sr');
hold on; loglog(normOut(:, 2), normOut(:, 4), '-ob'); 
hold on; loglog([2e-1, 4e-1], [1e-5, 16e-5], '-.r', 'LineWidth', 1.0); 
hold on; loglog([2e-1, 4e-1], [5e-9, 5*64e-9], '-.b', 'LineWidth', 1.0); 
legend('MIM4-COMP4', 'MIM6-COMP6', 'MIM4-RK4', ...
    'order = 4', 'order = 6', 'location', 'best') %
title('Numerical Convergence, 1D Boussinesq Eq')
xlabel('\Delta x'); ylabel('|| U - u_{exact} ||_{\infty}')
set(gca,'FontSize',10); grid on; 




function [uCOMP4, ECOMP4, tCOMP4] = COMP4(NElem,dh,dt,tEnd,u0, cV)

% 4th order Composition method
% Hairer et al, Geometric Integration, pp152, Sec V.3, eq. (3.6)

U = u0(1:NElem+2,1); Ut = u0(NElem+3:end,1); 
V = zeros(size(U));
uCOMP4(:,1) = [U; Ut]; % initial value
I1 = speye(NElem+2,NElem+2);
D = div(4, NElem, dh); G = grad(4, NElem, dh);
L = D*G; 

alp1 = (146 + 5*sqrt(19))/540;  bet5 = alp1;
alp2 = (-2 + 10*sqrt(19))/135;  bet4 = alp2; 
alp3 = 1/5;                     bet3 = alp3;
alp4 = (-23 - 20*sqrt(19))/270; bet2 = alp4;
alp5 = (14 - sqrt(19))/108;     bet1 = alp5; 

E1 = 0.5*(V'*V); E2 = 0.5*(Ut'*Ut);
E3 = 0.5*(U'*U); E4 = 1/3*(U.^2)'*U; 

E0 = dh*(E1 + E2 + E3 + E4);


ECOMP4(1,1) = E0; 
tCOMP4(1,1) = 0; % initial time    
iCount = 2; 
t = dt; 

V = Ut; 

    while t <= tEnd
        p1 = U + bet1*dt*V; 
        q1 = V + (I1 - L)\((bet1 + alp1)*dt*(L*(p1 + p1.^2)));         
        
        % Step 2        
        p2 = p1 + (bet2 + alp1)*dt*q1; 
        q2 = q1 + (I1 - L)\((bet2 + alp2)*dt*(L*(p2 + p2.^2))); 
    
    %
        % Step 3    
        p3 = p2 + (bet3 + alp2)*dt*q2; 
        q3 = q2 + (I1 - L)\((bet3 + alp3)*dt*(L*(p3 + p3.^2))); 
    
    
        % Step 4    
        p4 = p3 + (bet4 + alp3)*dt*q3; 
        q4 = q3 + (I1 - L)\((bet4 + alp4)*dt*(L*(p4 + p4.^2))); 

        % Step 5    
        p5 = p4 + (bet5 + alp4)*dt*q4; 
        q5 = q4 + (I1 - L)\((bet5 + alp5)*dt*(L*(p5 + p5.^2))); 
        
        p5 = p5 + alp5*dt*q5; 

        uNew = [p5; q5];                  
        
        Uin = uNew(1:NElem+2,1); Vin = uNew(NElem+3:end,1);
        Ut1 = cV*Vin;  
        E1 = 0.5*(Vin'*Vin); E2 = 0.5*(Ut1'*Ut1); 
        E3 = 0.5*(Uin'*Uin); E4 = 1/3*(Uin.^2)'*Uin; 

        ECOMP4(iCount,1) = dh*(E1 + E2 + E3 + E4) ; 


        uCOMP4(:,iCount) = uNew;   
        U = uCOMP4(1:NElem+2, iCount);
        V = uCOMP4(NElem+3:end, iCount);

        tCOMP4(iCount,1) = t;
        
        iCount = iCount+1; 
        t = t + dt;              
        
    end

end


function [uCOMP4, ECOMP4, tCOMP4] = COMP61(NElem,dh,dt,tEnd,u0, cV)

% 6th order Composition Method
% Hairer et al, Geometric Integration, pp152, Sec V.3
% Eq (3.9), (3.11)

U = u0(1:NElem+2,1); Ut = u0(NElem+3:end,1); 
V = zeros(size(U)); 
uCOMP4(:,1) = [U; V]; % initial value

I1 = speye(NElem+2,NElem+2);
D = div(6, NElem, dh); G = grad(6, NElem, dh);
L = D*G; 

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


E1 = 0.5*(V'*V); E2 = 0.5*(Ut'*Ut); E3 = 0.5*(U'*U);
E4 = 1/3*(U.^2)'*U; 
E0 = dh*(E1 + E2 + E3 +E4);
ECOMP4(1,1) = E0; 
ECOMP4(1,1) = 0; 
tCOMP4(1,1) = 0; % initial time    
iCount = 2; 
t = dt; 
V = Ut; 

    while t <= tEnd
        p1 = U + bet1*dt*V; 
        q1 = V + (I1 - L)\((bet1 + alp1)*dt*(L*(p1 + p1.^2)));         
        
        % Step 2        
        p2 = p1 + (bet2 + alp1)*dt*q1; 
        q2 = q1 + (I1 - L)\((bet2 + alp2)*dt*(L*(p2 + p2.^2))); 
    
    %
        % Step 3    
        p3 = p2 + (bet3 + alp2)*dt*q2; 
        q3 = q2 + (I1 - L)\((bet3 + alp3)*dt*(L*(p3 + p3.^2))); 
    
    
        % Step 4    
        p4 = p3 + (bet4 + alp3)*dt*q3; 
        q4 = q3 + (I1 - L)\((bet4 + alp4)*dt*(L*(p4 + p4.^2))); 

        % Step 5    
        p5 = p4 + (bet5 + alp4)*dt*q4; 
        q5 = q4 + (I1 - L)\((bet5 + alp5)*dt*(L*(p5 + p5.^2))); 

        % Step 6    
        p6 = p5 + (bet6 + alp5)*dt*q5; 
        q6 = q5 + (I1 - L)\((bet6 + alp6)*dt*(L*(p6 + p6.^2))); 

        % Step 7
        p7 = p6 + (bet7 + alp6)*dt*q6; 
        q7 = q6 + (I1 - L)\((bet7 + alp7)*dt*(L*(p7 + p7.^2))); 

        p7 = p7 + alp7*dt*q7; 

        uNew = [p7; q7];                  
        
        Uin = uNew(1:NElem+2,1); Vin = uNew(NElem+3:end,1);
        Utin = cV*Vin;   
        E1 = 0.5*(Vin'*Vin); E2 = 0.5*(Utin'*Utin); 
        E3 = 0.5*(Uin'*Uin); 
        E4 = 1/3*(Uin.^2)'*Uin; 

        ECOMP4(iCount,1) = dh*(E1 + E2 + E3 +E4); 

        uCOMP4(:,iCount) = uNew;   
        U = uCOMP4(1:NElem+2, iCount);
        V = uCOMP4(NElem+3:end, iCount);

        tCOMP4(iCount,1) = t;
        
        iCount = iCount+1; 
        t = t + dt;              
        
    end

end


function [uCOMP4, ECOMP4, tCOMP4] = COMP62(NElem,dh,dt,tEnd,u0, cV)

% 6th order Composition Method
% Hairer et al, Geometric Integration, pp152, Sec V.3
% Eq (3.9), (3.12)

U = u0(1:NElem+2,1); V = u0(NElem+3:end,1); 
uCOMP4(:,1) = [U; V]; % initial value

D = div(6, NElem, dh);
G = grad(6, NElem, dh);
L = D*G; 

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
E0 = 0.5*dh*(ut'*ut + V'*V + 2*(1. - cos(U))'*ones(size(U,1), 1));
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
        ECOMP4(iCount,1) = 0.5*dh*(ut1'*ut1 + q9'*q9 + 2*(1.- cos(p9))'*ones(size(p9,1), 1)) - E0; 
        uCOMP4(:,iCount) = uNew;   
        U = uCOMP4(1:NElem+2, iCount);
        V = uCOMP4(NElem+3:end, iCount);

        tCOMP4(iCount,1) = t;
        
        iCount = iCount+1; 
        t = t + dt;              
        
    end

end






function [uRRK, tRRK, gamRRK, eRRK] = RRK(RKFlag,NElem,dh,dt,tEnd,Uin0, cV)
    
% Relaxation RK4
     
    I1 = speye(NElem+2,NElem+2);
    D = div(4, NElem, dh); G = grad(4, NElem, dh);
    L = cV^2*D*G; 
                
    C(1,1) = 0;    C(2,1) = 1/2;
    C(3,1) = 1/2;    C(4,1) = 1;

    A(2,1) = 1/2;
    A(3,1) = 0;   A(3,2) = 1/2;
    A(4,1) = 0;   A(4,2) = 0;     A(4,3) = 1;     

    B(1,1) = 1/6;    B(2,1) = 1/3;
    B(3,1) = 1/3;    B(4,1) = 1/6;
     
    tRRK(1,1) = 0; % initial time    
    gamRRK(1,1) = 1;
    
    u0 = Uin0(1:NElem+2,1); ut0 = Uin0(NElem+3:end,1);     
    v0 = zeros(size(u0)); 

    E1 = 0.5*(v0'*v0); E2 = 0.5*(ut0'*ut0);
    E3 = 0.5*(u0'*u0); E4 = 1/3*(u0.^2)'*u0; 

    E0 = dh*( E1 + E2 + E3 + E4 );
    eRRK(1,1) = E0; 
    uRRK(:,1) = [u0; ut0]; % initial value
    y = [u0; ut0];
    
    iCount = 2; 
    t = dt; 
   
    while t <= tEnd

        z1 = y; 
        z1U = z1(1:NElem+2); z1V = z1(NElem+3:end); 
        [k1] = [z1V; 
                (I1 - L)\L*(z1U + z1U.^2)]; %
        
        % Step 2        
        z2 =  y + dt*A(2,1)*k1;       
        z2U = z2(1:NElem+2); z2V = z2(NElem+3:end); 
        [k2] = [z2V; 
                (I1 - L)\L*(z2U + z2U.^2)]; %

        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
        z3U = z3(1:NElem+2); z3V = z3(NElem+3:end); 
        [k3] = [z3V; 
                (I1 - L)\L*(z3U + z3U.^2)]; %      

        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        z4U = z4(1:NElem+2); z4V = z4(NElem+3:end); 
        [k4] = [z4V; 
                (I1 - L)\L*(z4U + z4U.^2)]; %

       
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
        
        Uin = uNew(1:NElem+2,1); Vin = uNew(NElem+3:end,1);
        ut1 = cV*Vin; 

        E1 = 0.5*(Vin'*Vin); E2 = 0.5*(ut1'*ut1);
        E3 = 0.5*(Uin'*Uin); E4 = 1/3*(Uin.^2)'*Uin;

        eRRK(iCount,1) = dh*(E1 + E2 + E3 + E4) ; 

        iCount = iCount+1; 
        t = t + gam*dt;           
        
    end   

end






function [u0, ut0] = u0Func(xGrid, CaseNo)

switch CaseNo
    case 1
    % case 1 - refer Li, J. Comp Phy, 2020
        al = 0.5; x0 = 0;
        bt = sqrt(1+2*al/3);
        Par = sqrt(al/6)*(xGrid - x0)/bt;
        u0 = al*sech(Par).^2; 
        Ex1 = 2*al*sqrt(al/6)*sech(Par).^2; 
        ut0 = Ex1.*tanh(Par); 
%         ut0 = zeros(size(u0,1),1); 

    
end

end


function [yR] = tInterp(tR, uR, tEnd, NElem)

    t1 = tR(end-1, 1); t2 = tR(end, 1);
    u1 = uR(1:NElem+2, end-1); u2 = uR(1:NElem+2, end); 
    yR = u1 + (u2 - u1)*(tEnd - t1)/(t2 - t1); 

end

function [uOut, utOut] = uExact(t, xGrid, CaseNo)

switch CaseNo    
    case 1  
        al = 0.5; x0 = 0;
        bt = sqrt(1+2*al/3);
        Par = sqrt(al/6)*(xGrid - bt*t - x0)/bt;
       uOut = al*sech(Par).^2; 
        Ex1 = 2*al*sqrt(al/6)*sech(Par).^2; 
        utOut = Ex1.*tanh(Par); 
        
       

end

          

end
