% The program solves the VFE for helical-circular-polygons in the Minkowski
% space
% Date: November 21, 2019
% -----------------------------------------

clear 

tic 
for kk = 9
M = 3 ; 
p1 = 1 ; 
q1 = 500*2^(kk-1); 
q1 = 3; 
p1 = 10327; q1 = 27993; 
N = 2^5 * M ; 
b = 1+1e-5 ; 
theta0 = 2*atan(b*tan(pi/M)) ; 
a = sqrt(-1+b^2) ; 
rho0 = 2 * asinh(a*sin(pi/M)) ;
c0 = sqrt(2*log(cosh(rho0/2))/pi) ; 
    
if p1 == 0 
    q = 1 ;
    p = 0 ; 
end

if ne(gcd(p1,q1),1) 
    p = p1 / gcd(p1,q1) ;
    q = q1 / gcd(p1,q1) ;
else 
    p = p1 ; 
    q = q1 ;
end

% ---------  file name 
Matfile = ['Alg_sol_VFE_circ_hel_M' num2str(M) '_b' num2str(b*1e1) 'p' num2str(p) '_q' num2str(q)] ; 

% --------------------------------------------------------
%     Initial values of angle rho and theta at time t = 0 
% --------------------------------------------------------   
 
%     if theta_0 >= 2*pi/M
%         Error('b>=1') 
%     end
    
    
    for m = 0 : q-1
        if mod(q,2) == 1
            t(m+1) = -1i*log(Gauss_sum(-p,m,q) / sqrt(q)) ; 
            rho(m+1) = 2*acosh(cosh(rho0/2)^(1/q)) ;
        elseif mod(q/2-m,2) == 0 
            t(m+1) = -1i*log(Gauss_sum(-p,m,q) / sqrt(2*q)) ; 
            rho(m+1) = 2*acosh(cosh(rho0/2)^(2/q))      ;
        else 
            t(m+1) = 0 ;
            rho(m+1) = 0;
        end
    end

    for k = 0 : M-1
        for m = 0 : q-1
            th(q*k+m+1) = theta0 *(k+ m/q) + t(m+1) ;
            rho_q(q*k+m+1) = rho(m+1) ; 
        end
    end
    

th = real(th) ;
th = th + p * theta0^2 / (2*pi*q) ; 

t = real(t) ; 
j = 1 ; 

% Computation of Rotation matrix 
% % Adding the addtional phase term 
% th = th+(theta_0^2/(2*pi))*(p/q); 

for j = 1 : M*q 
    R(:,:,j) = [cosh(rho_q(j))     sinh(rho_q(j))*cos(th(j))     sinh(rho_q(j))*sin(th(j)) ; 
        sinh(rho_q(j))*cos(th(j))    cosh(rho_q(j))*cos(th(j))^2+sin(th(j))^2      (cosh(rho_q(j))-1)*cos(th(j))*sin(th(j)) ; 
        sinh(rho_q(j))*sin(th(j))   (cosh(rho_q(j))-1)*cos(th(j))*sin(th(j))   cosh(rho_q(j))*sin(th(j))^2+cos(th(j))^2] ; 
end

% Product of M*q matrices 
Prd_R = eye(3) ; 
for j = 1 : M*q 
    Prd_R = R(:,:,j) * Prd_R ; 
end 

% Computation of Frenet frame {T,e1,e2} and tangent vectors 
Te1e2 = eye(3) ; 
T(1,:) = Te1e2(1,:) ; 

for k = 1 : M*q 
    Te1e2 = R(:,:,k) * Te1e2 ;
    T(k+1,:) = Te1e2(1,:) ; 
end 
       
T1 = T ; % Tangent vectors for plotting a closed curve 
T = T(1:end,:) ; 

% ---------------------------------------------
% Constructing X from the algorithm written in the draft
% It works for all kind of theta_0 or b values
% ---------------------------------------------
X(1,:) = [ 0 0 0 ] ; 

% Compute spq
spq = 2*theta0 * p / (M*q) ; 

X(2,:) = X(1,:) + spq * T(1,:) ; 
for j = 2 : (M*q)+1
    X(j+1,:) = X(j,:) + (2*pi/(M*q)) * T(j,:) ;
end


XX= X(2:end,:) ; 
TT = T(2:end,:) ; 
% ------------------------------------------------
% Rotation so that X satisfies X(2*pi)-X(0) = 2*pi*b 
% ------------------------------------------------

v1 = XX(end,:) - XX(1,:) ; % time-like vector
v1 = v1 / sqrt(v1(1)^2-v1(2)^2-v1(3)^2) ;
% Computation of angle with z axis and rotation axis

ang1 = acosh(v1(1)) ; % angle between two time-like vectors
u = cross(v1,[1 0 0]) ; u(1) = -u(1);  % axis
ax1 = u / sqrt(-u(1)^2+u(2)^2+u(3)^2) ; 

K1 = [0 ax1(3) -ax1(2) ; ax1(3) 0 -ax1(1) ; -ax1(2) ax1(1) 0 ] ; % Semi-skew symmetric 
% K1 = [0 ax1(3) ax1(2) ; ax1(3) 0 ax1(1) ; ax1(2) -ax1(1) 0 ] ; % Semi-skew symmetric 
R1 = eye(3) + sinh(ang1) * K1 + (cosh(ang1)-1) * K1^2 ; 
        
RX = (R1 * XX.').' ; 
RT = (R1 * TT.').' ; 

RT1= (R1 * T1.').' ; 

R2X =  RX ; 
R2T =  RT ; 

% Mean of vertices ONLY

R2X(:,1:2) = R2X(:,1:2) - mean(R2X(2:end-1,1:2));
% R2X(:,3) = R2X(:,3) - mean(R2X(2:end-1,3)); 
% R2X = R2X - mean(R2X);

c_M_alg = -2*log(cos(rho0/2))/(pi*tan(pi/M)/M) ; 

% R2X(:,3) = R2X(:,3) + c_M_alg*2*pi/M^2 ; 

% X3M(p1+1) = R2X(2,3) ; 

% Interpolation 
XXXX = zeros(N, 3);
TTTT = zeros(N, 3);

ssqM = spq + (2 * pi * (0:q*M) / (q * M) );

index1 = floor((q * M) * (0:N-1) / N) + 1;
index2 = floor((q * M) * (0:N-1) / N) + 2;

ss = spq + (2 * pi * (0:N-1) / N );
ss1 = ssqM(index1);
ss2 = ssqM(index2);
XXX = R2X ; TTT = R2T ; 
for m = 1:N
    XXXX(m, :) = ((ss(m) - ss2(m)) * XXX(index1(m), :) + (ss1(m) - ss(m)) * XXX(index2(m), :)) / (ss1(m) - ss2(m));
    TTTT(m, :) = TTT(index1(m), :);
end


save(Matfile, 'T','X','M', 'b','p','q', 'R2X','R2T','-v7.3' ) ;

%  h0 = mean((0:2*pi/M:2*pi)*b); % mean(X3) at t=0
% h_tpq = h0 + c_M_alg*2*pi/M^2*p/q ; % mean(X3) at t = tpq 

end 
return ;

 h0 = pi*b*(M-1)/M ; % mean(X3) at t=0 for M polygon
 h_tpq = mean(R2X(2:q+1,3));
 
 diff = R2T(end,:) - R2T(2,:);
 ds = 2*pi/(M*q);
t_1q = 2*pi/(M^2*q) ;
c0_num = sqrt(t_1q) * norm(diff) / (2*ds);
c0;
abs(c0-c0_num)












