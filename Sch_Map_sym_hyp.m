% In this program we solve VFE and Schrodinger Map equation for a regular
% helical polygon in the Hyperbolic plane. We make use of the Spatial symmetry of T(s,t) and 
% reduce the compuation cost by reducing the number of elements to N/M
% Hyperbolic metric : d(x,y) = -x1y1+x2y2+x3y3
% Date: June 11, 2018 

close all;
clear  

% Sch_Map

% function Sch_Map
tic

	% function call
% 	cleanupObj = onCleanup(@cleanMeUp);

    M = 3 ; 
    n = -2 ; 
    q = 2*3^2*5 ; p = 1 ;  
%     q = 2^4*3^2*5*7 ; p = 1 ;  
    N = 512* 2^n * M; 
%     N = M * q ; 
%     N = 2^(3+n) * 3 * 5 * M  ; 
    N_t = 151200 * 4^n ;
%         N_t = 151200 * 4^n * 1.5   ; 
%     N_t = 4^(n+3)  * 3^3 * 5 ; 
    l = 1 ;                 % # of loops
%     theta = pi/2 ; 
%     b = tan(theta/2)/tan(pi/M) ;             % For hyperbolic case
    b = 1.2 ; 
    L = 2*pi ; 
    s =  (0 : 1/N : 1-1/N ) * L ; s = s.' ; 
    h = s(2) - s(1) ; 
    
   T_run =  2 * pi / M^2 ; 
   
   dt = T_run / N_t ;
    
    k1 = 2i * pi * [0 : N/2-1 -N/2 : -1 ] / L ; k1 = k1.' ; 
    
% k1 = 2i * pi * [0:N/2-1 -N/2:-1] / L;
     
% Initial data 
    z0 = zeros(N,1) ; 
    for j = 0 : M-1
        z0( j* N/M + (1:N/M) ) = exp(2i*pi* j / M) ; 
    end 
    j = 0 : 1 ; 
    x_0 =  (pi * exp(1i * pi * (2 * j - 1 ) / M )) ./ (M * sin(pi / M) )   ;  % Initial data
     
    X00 = linspace(x_0(1), x_0(2), N/M + 1);
    Z0 = zeros(N,1);

    Z0(1:N/M) = X00(1:N/M);
    for m = 1: M-1
        Z0(m * N/M + (1:N/M)) = exp(2i * pi * m / M) * Z0(1:N/M); 
    end

    T10 = b * ones(N,1) ; T20 = sqrt(b^2-1) * imag(z0) ; T30 = sqrt(b^2-1) * real(z0) ; 
    X10 = b * s ; X20 = -sqrt(b^2-1) * real(Z0); X30 = sqrt(b^2-1) * imag(Z0); 
    % Storing the initial data
    Xfull = zeros(N/M,q+1) ; XX1full = zeros(N/M,q+1) ; 
    Tfull = zeros(N/M,q+1) ; TT1full = zeros(N/M,q+1) ; 
    
    T1 = T10; T2 = T20; T3 = T30 ; 
    X1 = X10 ; X2 = X20 ; X3 = X30 ; 
    
    T = T2 + 1i*T3 ; X = X2 + 1i*X3 ; 
    
    T = T(1:N/M) ; X = X(1:N/M) ; 
    T1 = T1(1:N/M) ; X1 = X1(1:N/M) ; 
    
    Xfull(:,1) = X ; XX1full(:,1) = X1 ; 
    Tfull(:,1) = T; TT1full(:,1) = T1 ; 
    z1(1,:) = [X1(1) X2(1) X3(1)] ;
        
    % X_mean 

    X_mean = zeros(N_t,1) ; X1_mean = zeros(N_t,1) ; 

%   X_mean(1) = h * (0.5 * X(1) + sum(X(2:end-1)) + 0.5 * X(end)) / (L/M)  ;
%   X1_mean(1) = h * (0.5 * X1(1) + sum(X1(2:end-1)) + 0.5 * X1(end)) / (L/M)  ;
    X_mean(1) = sum(X) / (N/M) ; 
    X1_mean(1) = sum(X1) / (N/M) ; 
    
% Stability condition     
    if dt/h^2 > 0.29
        error('Step size is too large') 
    end
    
%   Frequencies for N/M elements 

    kk1 = 2* pi * [0 : (N)/(2*M)-1 -(N)/(2*M) : -1 ] / L ; kk1 = kk1.' ;
    k2 = M*kk1+1 ; 
    k3 = M* kk1 ; 

% Rotation R 
    j = 0 : N/M-1 ;
    R = exp(-2i*pi*j/N) ; R = R.' ; 

%     plt = animatedline ; 

for r = 1 : N_t
    
    if isnan(T1(1))
        r
        r * dt
        error('Error!');
    end
    
    T1_ =  fft(T1) ; T1_(abs(T1_) < eps) = 0 ; T1s = real(ifft( T1_ * 1i .* k3 )) ; T1ss = real(ifft( T1_ .* (1i * k3).^2 )) ; 
    Tr = T .* R ; Tr_ =  fft(Tr) ; Tr_(abs(Tr_)/N < eps) = 0 ; Ts = ifft( Tr_ * 1i .* k2 ) .* conj(R) ; Tss = ifft( Tr_ .* (1i * k2 ).^2 ) .* conj(R) ; 
    
    AT1 = -(real(T) .* imag(Tss) - imag(T) .* real(Tss)) ; 
    AT = 1i * (-T1ss .* T  + T1 .* Tss) ; 
    AX1 = -(real(T) .* imag(Ts) - imag(T) .* real(Ts)) ; 
    AX = 1i * (-T1s .* T  + T1 .* Ts) ; 
    T1aux = T1 + 0.5 * dt * AT1 ; 
    Taux = T + 0.5 * dt * AT ; 
    
    T1_ =  fft(T1aux) ; T1_(abs(T1_) < eps) = 0 ; T1s = real(ifft( T1_ * 1i .* k3 )) ; T1ss = real(ifft( T1_ .* (1i * k3).^2 )) ; 
    Tr = Taux .* R ; Tr_ = fft(Tr) ; Tr_(abs(Tr_)/N < eps) = 0 ; Ts = ifft( Tr_ * 1i .* k2 ) .* conj(R) ; Tss = ifft( Tr_ .* (1i * k2 ).^2 ) .* conj(R) ; 
    
    BT1 = -(real(Taux) .* imag(Tss) - imag(Taux) .* real(Tss)) ; 
    BT = 1i * (-T1ss .* Taux  + T1aux .* Tss) ; 
    BX1 = -(real(Taux) .* imag(Ts) - imag(Taux) .* real(Ts)) ; 
    BX = 1i * (-T1s .* Taux  + T1aux .* Ts) ; 
    T1aux = T1 + 0.5 * dt * BT1 ; 
    Taux = T + 0.5 * dt * BT ; 
    
    
    T1_ = fft(T1aux) ; T1_(abs(T1_) < eps) = 0 ; T1s = real(ifft( T1_ * 1i .* k3 )) ; T1ss = real(ifft( T1_ .* (1i * k3).^2 )) ; 
    Tr = Taux .* R ; Tr_ = fft(Tr) ; Tr_(abs(Tr_)/N < eps) = 0 ; Ts = ifft( Tr_ * 1i .* k2 ) .* conj(R) ; Tss = ifft( Tr_ .* (1i * k2 ).^2 ) .* conj(R) ; 
    
    CT1 = -(real(Taux) .* imag(Tss) - imag(Taux) .* real(Tss)) ; 
    CT = 1i * (-T1ss .* Taux  + T1aux .* Tss) ; 
    CX1 = -(real(Taux) .* imag(Ts) - imag(Taux) .* real(Ts)) ; 
    CX = 1i * (-T1s .* Taux  + T1aux .* Ts) ; 
    Taux = T +  dt * CT ; 
    T1aux = T1 +  dt * CT1 ; 
    
    T1_ = fft(T1aux) ; T1_(abs(T1_) < eps) = 0 ; T1s = real(ifft( T1_ * 1i .* k3 )) ; T1ss = real(ifft( T1_ .* (1i * k3).^2 )) ; 
    Tr = Taux .* R ; Tr_ = fft(Tr) ; Tr_(abs(Tr_)/N < eps) = 0 ; Ts = ifft( Tr_ * 1i .* k2 ) .* conj(R) ; Tss = ifft( Tr_ .* (1i * k2).^2 ) .* conj(R) ; 
    
    DT1 = -(real(Taux) .* imag(Tss) - imag(Taux) .* real(Tss)) ; 
    DT = 1i * (-T1ss .* Taux  + T1aux .* Tss) ; 
    DX1 = -(real(Taux) .* imag(Ts) - imag(Taux) .* real(Ts)) ; 
    DX = 1i * (-T1s .* Taux  + T1aux .* Ts) ; 
    
    T1 = T1 + dt * (AT1 + 2 * BT1 + 2 * CT1 + DT1) / 6 ; 
    T = T + dt * (AT + 2 * BT + 2 * CT + DT) / 6 ; 
    X1 = X1 + dt * (AX1 + 2 * BX1 + 2 * CX1 + DX1) / 6 ; 
    X = X + dt * (AX + 2 * BX + 2 * CX + DX) / 6 ; 
    
    
%     Tnorm = sqrt(real(T).^2+imag(T).^2 + T3.^2) ; 
%     T = T ./ Tnorm ; T3 = T3 ./ Tnorm ; 

    Tnorm = sqrt( T1.^2 -real(T).^2-imag(T).^2 ) ; 

    T = T ./ Tnorm ; T1 = T1 ./ Tnorm ; 

    if mod(r,N_t*p/q) == 0 
        Xfull(:,p+1) = X ;  XX1full(:,p+1) = X1 ; 
        Tfull(:,p+1) = T ;  TT1full(:,p+1) = T1 ; 
        p = p + 1 ; 
%         plot3(X1,X2,X3,'o-'), drawnow 
%         plot3(T1,real(T),imag(T)) , drawnow
%         plot(T2./(1+T1),T3./(1+T1)), drawnow 
    end 
%     figure, hold on, for r = 0:.01:100, plot3(cos(r), sin(r), r), drawnow, end

    z1(r+1,:) = [X1(1) real(X(1)) imag(X(1)) ] ;
% 
%     X_mean(r+1) = h * (0.5 * X(1) + sum(X(2:end-1)) + 0.5 * X(end)) / (L/M)  ;
%     X1_mean(r+1) = h * (0.5 * X1(1) + sum(X1(2:end-1)) + 0.5 * X1(end)) / (L/M)  ;
    
    X1_mean(r+1) = sum(X1) / (N/M) ;
    X_mean(r+1) = sum(X) / (N/M) ;
end
    
T1full(:,:) = zeros(l*N,q+1) ; X1full = zeros(l*N,q+1); 
TTfull = zeros(l*N,q+1) ; XXfull = zeros(l*N,q+1); 
    for k = 1 : q+1
    
        for m = 0 : l*M-1
            TTfull(m * N/M + (1:N/M),k) = exp(2i * pi * m / M) * Tfull(:,k); 
            T1full(m * N/M + (1:N/M),k) =  TT1full(:,k); 
   
            XXfull(m * N/M + (1:N/M),k) = exp(2i * pi * m / M) * Xfull(:,k); 
            X1full(m * N/M + (1:N/M),k) = 2*pi*b*m/M + XX1full(:,k); 
        end
    
    end
    T2full = real(TTfull) ; T3full = imag(TTfull) ;
    X2full = real(XXfull) ; X3full = imag(XXfull) ;




% Calculation for h and c_M
% h_T = M * sum(X3) ; 
c_M = (X1_mean(end)-X1_mean(1)) / T_run ; 


save Sch_Map_fft_hyp_M6_2N7_dt_hyp_theta_pi2.mat X1full X2full X3full T1full T2full T3full X_mean X1_mean z1 L  M dt b c_M


toc
return ; 

% % fires when main function terminates
%     function cleanMeUp()
%         % saves data to file (or could save to workspace)
%         fprintf('saving on cleanup...\n');
% %         Sch_Map_fft_sym_M3_NMq_dt_euc_b_0_q24_32_5.mat X1full X2full X3full T1full T2full T3full X_mean X3_mean z1 L M dt b c_M  
% %         filename = [datestr(now,'yyyy-mm-dd_HHMMSS') '.mat'];
%           save Sch_Map_fft_sym_M3_NMq_dt_euc_b_0_q24_32_5.mat  X_mean X3_mean z1 L M dt b c_M  
%     end
%     
% end 


