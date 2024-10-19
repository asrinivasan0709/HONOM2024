clc; clear; close all;
% Symplectic Scheme, 1D Wave eqn, Composition method
% ASrinivasan, 1Apr2024
% Uses fourth and sixth order mimetic operators

addpath('...\mole_MATLAB')
Mult = 2; 
m = Mult*100; % xmin = -40; xmax = 40; 
xmin = -15; xmax = 15; 
dx = (xmax-xmin)/m; tEnd = 12; cV = 1; 
CFL = 0.5; dt = CFL*dx/cV; 

% initial condition
xGrid = [xmin xmin+dx/2:dx:xmax-dx/2 xmax]';
xNod = [xmin:dx:xmax]';
CaseNo = 4;  
[u0, v0] = u0Func(xGrid, CaseNo); u0 = [u0; v0]; 

% Energy Calcs
[uCOMP4, ECOMP4, tCOMP4] = COMP4(m,dx,dt,tEnd,u0,cV); 
[uCOMP6, ECOMP6, tCOMP6] = COMP61(m,dx,dt,tEnd,u0,cV); 
[uRRK, tRRK, gamRRK, eRRK] = RRK(0,m,dx,dt,tEnd,u0,cV); 
uEx = uExact(tRRK(end), xGrid, CaseNo); 

normCalc = [norm(uRRK(1:m+2, end) - uEx, 2), ...
                    norm(uCOMP4(1:m+2, end) - uEx, 2), ...
                    norm(uCOMP6(1:m+2, end) - uEx, 2)]; 


figure; 
plot(xGrid, uCOMP4(1:m+2, 1), '--k', 'LineWidth', 1);
hold on; plot(xGrid, uCOMP4(1:m+2, end), '-.k', 'LineWidth', 2);
hold on; plot(xGrid, uCOMP6(1:m+2, end), 'sr', 'LineWidth', 2);
hold on; plot(xGrid, uRRK(1:m+2, end), '--b', 'LineWidth', 1.2);
hold on; plot(xGrid, uEx, '-.r', 'LineWidth', 1.5);
title('u vs x, 1D Wave Eq'); 
legend('initial', 'MIM4-COMP4', 'MIM6-COMP6', 'MIM4-RK4', 'exact', 'location', 'best'); 
xlabel('x'); ylabel('u(t,x)'); set(gca,'FontSize',10)

figure; 
plot(xGrid, uCOMP4(1:m+2, 1), '--k', 'LineWidth', 1);
hold on; plot(xGrid, uCOMP4(1:m+2, end), '-.k', 'LineWidth', 2);
hold on; plot(xGrid, uCOMP6(1:m+2, end), 'sr', 'LineWidth', 2);
hold on; plot(xGrid, uRRK(1:m+2, end), '--b', 'LineWidth', 1.2);
hold on; plot(xGrid, uEx, '-.r', 'LineWidth', 1.5);
title('u vs x, 1D Wave Eq'); 
legend('initial', 'MIM4-COMP4', 'MIM6-COMP6', 'MIM4-RK4', 'exact', 'location', 'best'); 
xlabel('x'); ylabel('u(t,x)'); set(gca,'FontSize',10); xlim([9.5, 13.5]); ylim([-0.2 1.0]) ;

pn1 = 1; pn2 = 10; figure; 
plot(tCOMP4(pn1:pn2:end), abs(ECOMP4(pn1:pn2:end)/ECOMP4(1)), '-.k', 'LineWidth', 1.5); 
hold on; plot(tCOMP6(pn1:pn2:end), abs(ECOMP6(pn1:pn2:end)/ECOMP6(1)), 'sr', 'LineWidth', 1); 
hold on; plot(tRRK(pn1:pn2:end), abs(eRRK(pn1:pn2:end)/eRRK(1)), '--b', 'LineWidth', 1.2); 
title('Energy En/E0 vs time, 1D Wave Eq'); % set(gca, 'yscale', 'log');
ylim([0.9996 1.0001]); 
xlabel('time(s)'); ylabel('En/E0'); legend('MIM4-COMP4', 'MIM6-COMP6', 'MIM4-RK4', 'location', 'best')
set(gca,'FontSize',10); grid on; 

%%%% Convergence Calcs

VecM = [1,2,4,8,16,32]';

for i = 1:size(VecM,1)
    Mult = VecM(i); 
    m = Mult*64; xmin = -15; xmax = 15; 
    dx = (xmax-xmin)/m; tEnd = 1; cV = 1; 
    CFL = 1; dt = CFL*dx/cV; 

    % initial condition
    xGrid = [xmin xmin+dx/2:dx:xmax-dx/2 xmax]';
    xNod = [xmin:dx:xmax]';
    [u0, v0] = u0Func(xGrid, CaseNo); 
    u0 = [u0; v0]; 
    
    % Energy Calcs
    [uCOMP4, ECOMP4, tCOMP4] = COMP4(m,dx,dt,tEnd,u0,cV); 
    [uCOMP6, ECOMP6, tCOMP6] = COMP61(m,dx,dt,tEnd,u0,cV); 
    [uRRK, tRRK, gamRRK, eRRK] = RRK(0,m,dx,dt,tEnd,u0,cV); 
    uEx = uExact(tRRK(end), xGrid, CaseNo); 

    normCalc = [norm(uRRK(1:m+2, end) - uEx, 'inf'), ...
                    norm(uCOMP4(1:m+2, end) - uEx, 'inf'), ...
                    norm(uCOMP6(1:m+2, end) - uEx, 'inf')]; 

    normOut(i, :) = [m, dx, dt, normCalc];

end

% Calculate order of convergence
for i = 1:size(VecM, 1) - 1
    CvRK4(i, 1) = 1/log(2)*log(normOut(i,4)/normOut(i+1,4)); 
    CvCOMP4(i, 1) = 1/log(2)*log(normOut(i,5)/normOut(i+1,5)); 
    CvCOMP6(i, 1) = 1/log(2)*log(normOut(i,6)/normOut(i+1,6)); 
    
end

Cv = [[0;CvRK4], [0;CvCOMP4], [0;CvCOMP6]]; 

normTab = [normOut(:, 1:3), normOut(:,4), Cv(:,1), ...
            normOut(:,5), Cv(:,2), ...
            normOut(:,6), Cv(:,3)]; 
sympref('FloatingPointOutput',1);
LatTab = latex(sym(normTab))

function [uCOMP4, ECOMP4, tCOMP4] = COMP4(NElem,dh,dt,tEnd,u0, cV)

% 4th order Composition method
% Hairer et al, Geometric Integration, pp152, Sec V.3, eq. (3.6)

U = u0(1:NElem+2,1); V = u0(NElem+3:end,1); 
uCOMP4(:,1) = [U; V]; % initial value

D = div(4, NElem, dh);     G = grad(4, NElem, dh);
ID = interpDMat(4, NElem); IG = interpGMat(4, NElem);
D = D*ID; G = IG*G;

alp1 = (146 + 5*sqrt(19))/540;  bet5 = alp1;
alp2 = (-2 + 10*sqrt(19))/135;  bet4 = alp2; 
alp3 = 1/5;                     bet3 = alp3;
alp4 = (-23 - 20*sqrt(19))/270; bet2 = alp4;
alp5 = (14 - sqrt(19))/108;     bet1 = alp5; 

E0 = 0.5*dh*(U'*U + V'*V );
ECOMP4(1,1) = E0; 
tCOMP4(1,1) = 0; % initial time    
iCount = 2; 
t = dt; 

    while t <= tEnd
        p1 = U + bet1*dt*D*V; 
        q1 = V + (bet1 + alp1)*dt*(G*p1);    
        
        % Step 2        
        p2 = p1 + (bet2 + alp1)*dt*D*q1; 
        q2 = q1 + (bet2 + alp2)*dt*(G*p2); 
    
    %
        % Step 3    
        p3 = p2 + (bet3 + alp2)*dt*D*q2; 
        q3 = q2 + (bet3 + alp3)*dt*(G*p3); 
    
        % Step 4    
        p4 = p3 + (bet4 + alp3)*dt*D*q3; 
        q4 = q3 + (bet4 + alp4)*dt*(G*p4); 

        % Step 5    
        p5 = p4 + (bet5 + alp4)*dt*D*q4; 
        q5 = q4 + (bet5 + alp5)*dt*(G*p5);
    
        p5 = p5 + alp5*dt*D*q5; 

        uNew = [p5; q5];                  
        
        ECOMP4(iCount,1) = 0.5*dh*(p5'*p5 + q5'*q5 ); % - E0; 
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

D = div(6, NElem, dh);     G = grad(6, NElem, dh);
ID = interpDMat(6, NElem); IG = interpGMat(6, NElem);
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


E0 = 0.5*dh*(U'*U + V'*V );
ECOMP4(1,1) = E0; 
tCOMP4(1,1) = 0; % initial time    
iCount = 2; 
t = dt; 

    while t <= tEnd
        p1 = U + bet1*dt*D*V; 
        q1 = V + (bet1 + alp1)*dt*(G*p1);         
        
        % Step 2        
        p2 = p1 + (bet2 + alp1)*dt*D*q1; 
        q2 = q1 + (bet2 + alp2)*dt*(G*p2 ); 
    
    %
        % Step 3    
        p3 = p2 + (bet3 + alp2)*dt*D*q2; 
        q3 = q2 + (bet3 + alp3)*dt*(G*p3); 
    
    
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
        
        ECOMP4(iCount,1) = 0.5*dh*(p7'*p7 + q7'*q7 ); % - E0; 
        uCOMP4(:,iCount) = uNew;   
        U = uCOMP4(1:NElem+2, iCount);
        V = uCOMP4(NElem+3:end, iCount);

        tCOMP4(iCount,1) = t;
        
        iCount = iCount+1; 
        t = t + dt;              
        
    end

end


function [uRRK, tRRK, gamRRK, eRRK] = RRK(RKFlag,NElem,dh,dt,tEnd,u0, cV)
    
% Runge Kutta     

    D = div(4, NElem, dh);     G = grad(4, NElem, dh);
    ID = interpDMat(4, NElem); IG = interpGMat(4, NElem);
    D = D*ID; G = IG*G;
               
    C(1,1) = 0;    C(2,1) = 1/2;
    C(3,1) = 1/2;    C(4,1) = 1;

    A(2,1) = 1/2;
    A(3,1) = 0;   A(3,2) = 1/2;
    A(4,1) = 0;   A(4,2) = 0;     A(4,3) = 1;     

    B(1,1) = 1/6;     B(2,1) = 1/3;
    B(3,1) = 1/3;     B(4,1) = 1/6;
     
    uRRK(:,1) = u0; % initial value
    y = u0; 
    tRRK(1,1) = 0; % initial time    
    gamRRK(1,1) = 1;
    
    Uin0 = y(1:NElem+2,1); Vin0 = y(NElem+3:end,1); 
    E0 = 0.5*dh*(Uin0'*Uin0 + ...
            Vin0'*Vin0 ); 
    eRRK(1,1) = E0; % - E0; 
    
    iCount = 2; 
    t = dt; 
   
    while t <= tEnd

        z1 = y; 
        z1U = z1(1:NElem+2); z1V = z1(NElem+3:end); 
        [k1] = [D*z1V; 
                G*z1U ]; %

        % Step 2        
        z2 =  y + dt*A(2,1)*k1;       
        z2U = z2(1:NElem+2); z2V = z2(NElem+3:end); 
        [k2] = [D*z2V; 
                G*z2U ]; %


        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
        z3U = z3(1:NElem+2); z3V = z3(NElem+3:end); 
        [k3] = [D*z3V; 
                G*z3U ]; %  


        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        z4U = z4(1:NElem+2); z4V = z4(NElem+3:end); 
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
        
        Uin = uNew(1:NElem+2,1); Vin = uNew(NElem+3:end,1);
        eRRK(iCount,1) = 0.5*dh*(Uin'*Uin + ...
            Vin'*Vin ); % - E0;  
        
        iCount = iCount+1; 
        t = t + gam*dt;           
        
    end   

end



function [u0, v0] = u0Func(xGrid, CaseNo)

switch CaseNo
    case 1
    % case 1 - refer Mohebbi, Dehghan
        u0 = 4*sech(xGrid); 
        v0 = zeros(size(xGrid, 1), 1); 

    case 2
        ep = 0.; 
        u0 =  1/pi*sin(pi*(xGrid - ep)); 
        v0 = zeros(size(xGrid, 1), 1);  

    case 3
        % impulse, ref MA Sanchez et al
        for i = 1:size(xGrid, 1)
            bb = 0.; dl = 1;             
            P = (xGrid(i) - bb)/dl; 
            if abs(P) <= 0.5
                u0(i,1) = (2*P - 1)^10*(2*P + 1)^10; 
            else
                u0(i,1) = 0;
            end
        
        end

        v0 = -u0;  


    case 4
        A = 1; sigma = 0.5; mu = 0;         
        u0 = A*exp(-1/(sigma^2)*((xGrid-mu/2).^2));   
        v0 = -u0; 
    
end

end



function uOut = uExact(t, xGrid, CaseNo)

switch CaseNo
    case 1
        uOut = 4*atan(t*sech(xGrid)); 
        
    case 2
        ep = 0.; 
        uOut = 1/pi*sin(pi*(xGrid - ep)).*cos(pi*t); 

    case 3
        % impulse, case 4.4 of MA Sanchez et al
        for i = 1:size(xGrid, 1)
            bb = 0.; dl = 1;             
            P = (xGrid(i) - bb - t)/dl; 
            if abs(P) <= 0.5
                uOut(i,1) = (2*P - 1)^10*(2*P + 1)^10; 
            else
                uOut(i,1) = 0;
            end

        end

    case 4
        A = 1; sigma = 0.5; mu = 0;         
        uOut = A*exp(-1/(sigma^2)*((xGrid-mu/2 - t).^2));   
                  
end
         

end

