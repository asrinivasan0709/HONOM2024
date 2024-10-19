clc; clear; close all;
% Symplectic Scheme, 1D Sine Gordon eqn, Composition method
% ASrinivasan, 1Apr2024, Convergence

addpath('C:\Anand\Acer_Data\SDSU\MOLE\mole-master\mole_MATLAB')

VecM = [1,2,4,8,16]';

for i = 1:size(VecM,1)
    Mult = VecM(i); 
    m = Mult*64; xmin = -50; xmax = 50; 
    dx = (xmax-xmin)/m; tEnd = 15; cV = 1; 
    CFL = 0.25; dt = CFL*dx/cV; 
    
    normCalc = normFunc(m,dx,dt,tEnd,cV,xmin,xmax); 
    normOut(i, :) = [m, dx, dt, normCalc];

end

% Calculate order of convergence
for i = 1:size(VecM, 1) - 1
    CvRK4(i, 1) = 1/log(2)*log(normOut(i,4)/normOut(i+1,4)); 
    CvCOMP4(i, 1) = 1/log(2)*log(normOut(i,5)/normOut(i+1,5)); 
    CvCOMP6(i, 1) = 1/log(2)*log(normOut(i,6)/normOut(i+1,6)); 
    
end

Cv = [[0; CvRK4], [0; CvCOMP4], [0; CvCOMP6]];
normOut = [normOut, Cv]; 

figure;
loglog(normOut(:, 2), normOut(:, 5), '-*k');
hold on; loglog(normOut(1:4, 2), normOut(1:4, 6), '-sr');
hold on; loglog(normOut(:, 2), normOut(:, 4), '-ob'); 
hold on; loglog([2e-1, 4e-1], [5e-3, 5*16e-3], '-.r', 'LineWidth', 1.0); 
hold on; loglog([2e-1, 4e-1], [5e-6, 5*64e-6], '-.b', 'LineWidth', 1.0); 
legend('MIM4-COMP4', 'MIM6-COMP6', 'MIM4-RK4', ...
    'order = 4', 'order = 6', 'location', 'best') %
title('Numerical Convergence, 1D Sine Gordon Eq')
xlabel('\Delta x'); ylabel('|| U - u_{exact} ||_{\infty}')
set(gca,'FontSize',10); grid on; 



function normCalc = normFunc(m,dx,dt,tEnd,cV,xmin,xmax)

    % initial condition
    xGrid = [xmin xmin+dx/2:dx:xmax-dx/2 xmax]';
    xNod = [xmin:dx:xmax]';
    CaseNo = 2; % 2 works
    [u0, v0] = u0Func(xGrid, CaseNo); 
    u0 = [u0; v0]; 
    
    [uCOMP4, ECOMP4, tCOMP4] = COMP4(m,dx,dt,tEnd,u0,cV); 
    [uCOMP6, ECOMP6, tCOMP6] = COMP61(m,dx,dt,tEnd,u0,cV); 
    [uRRK, tRRK, gamRRK, eRRK] = RRK(0,m,dx,dt,tEnd,u0,cV); 
    uEx = uExact(tRRK(end), xGrid, CaseNo); 
    
    normCalc = [norm(uRRK(1:m+2, end) - uEx, 'inf'), ...
                    norm(uCOMP4(1:m+2, end) - uEx, 'inf'), ...
                    norm(uCOMP6(1:m+2, end) - uEx, 'inf')]; 

end





function [uCOMP4, ECOMP4, tCOMP4] = COMP4(NElem,dh,dt,tEnd,u0, cV)

% 4th order Composition method
% Hairer et al, Geometric Integration, pp152, Sec V.3, eq. (3.6)

U = u0(1:NElem+2,1); V = u0(NElem+3:end,1); 
uCOMP4(:,1) = [U; V]; % initial value

D = div(4, NElem, dh);
G = grad(4, NElem, dh);
L = D*G; 

alp1 = (146 + 5*sqrt(19))/540;  bet5 = alp1;
alp2 = (-2 + 10*sqrt(19))/135;  bet4 = alp2; 
alp3 = 1/5;                     bet3 = alp3;
alp4 = (-23 - 20*sqrt(19))/270; bet2 = alp4;
alp5 = (14 - sqrt(19))/108;     bet1 = alp5; 

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
        
        p5 = p5 + alp5*dt*q5; 

        uNew = [p5; q5];                  
        
        ut1 = cV*G*p5; 
        ECOMP4(iCount,1) = 0.5*dh*(ut1'*ut1 + q5'*q5 + 2*(1.- cos(p5))'*ones(size(p5,1), 1)) - E0; 
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

U = u0(1:NElem+2,1); V = u0(NElem+3:end,1); 
uCOMP4(:,1) = [U; V]; % initial value

D = div(6, NElem, dh);
G = grad(6, NElem, dh);
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

        p7 = p7 + alp7*dt*q7; 

        uNew = [p7; q7];                  
        
        ut1 = cV*G*p7; 
        ECOMP4(iCount,1) = 0.5*dh*(ut1'*ut1 + q7'*q7 + 2*(1.- cos(p7))'*ones(size(p7,1), 1)) - E0; 
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

% D = MimDMat(NElem, 6); 
% D = [zeros(1, NElem+1); D; zeros(1, NElem+1)]; 
% D = 1/dh*D; 
% G = MimGMat(NElem, 6); G = 1/dh*G; 
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








function [uRRK, tRRK, gamRRK, eRRK] = RRK(RKFlag,NElem,dh,dt,tEnd,u0, cV)
    
% Relaxation RK4
     
    D = div(4, NElem, dh);
    G = grad(4, NElem, dh);
    L = cV^2*D*G; 
    
               
    C(1,1) = 0;
    C(2,1) = 1/2;
    C(3,1) = 1/2;
    C(4,1) = 1;

    A(2,1) = 1/2;
    A(3,1) = 0;   A(3,2) = 1/2;
    A(4,1) = 0;   A(4,2) = 0;     A(4,3) = 1;     

    B(1,1) = 1/6; 
    B(2,1) = 1/3;
    B(3,1) = 1/3;
    B(4,1) = 1/6;

     
    uRRK(:,1) = u0; % initial value
    y = u0; 
    tRRK(1,1) = 0; % initial time    
    gamRRK(1,1) = 1;
    
    Uin0 = y(1:NElem+2,1); Vin0 = y(NElem+3:end,1); 
    ut = cV*G*Uin0;
    E0 = 0.5*dh*(ut'*ut + ...
            Vin0'*Vin0 + 2*(1. - cos(Uin0))'*ones(size(Uin0,1),1)); 
    eRRK(1,1) = 0; 
    
    iCount = 2; 
    t = dt; 
   
    while t <= tEnd

        z1 = y; 
        z1U = z1(1:NElem+2); z1V = z1(NElem+3:end); 
        [k1] = [speye(NElem+2,NElem+2)*z1V; 
                L*z1U - sin(z1U)]; %
        
        % Step 2        
        z2 =  y + dt*A(2,1)*k1;       
        z2U = z2(1:NElem+2); z2V = z2(NElem+3:end); 
        [k2] = [speye(NElem+2,NElem+2)*z2V; 
                L*z2U - sin(z2U)]; %
        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
        z3U = z3(1:NElem+2); z3V = z3(NElem+3:end); 
        [k3] = [speye(NElem+2,NElem+2)*z3V; 
                L*z3U - sin(z3U)]; %        
        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        z4U = z4(1:NElem+2); z4V = z4(NElem+3:end); 
        [k4] = [speye(NElem+2,NElem+2)*z4V; 
                L*z4U - sin(z4U)]; %   

       
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
        ut1 = cV*G*Uin;
       eRRK(iCount,1) = 0.5*dh*(ut1'*ut1 + ...
            Vin'*Vin + 2*(1. - cos(Uin))'*ones(size(Uin, 1),1)) - E0;  
        
        iCount = iCount+1; 
        t = t + gam*dt;           
        
    end   

end






function [u0, v0] = u0Func(xGrid, CaseNo)

switch CaseNo
    case 1
    % case 1 - refer Mohebbi, Dehghan
        u0 = zeros(size(xGrid, 1), 1); 
        v0 = 4*sech(xGrid); 

    case 2
        c = 0.5; eta = xGrid./sqrt(1-c^2);  
        u0 = 4*atan(exp(eta)); 
        Nr = 4*c*exp(eta); 
        Dr = sqrt(1-c^2)*(1 + exp(2*eta));
        v0 = -Nr./Dr; 
end

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
