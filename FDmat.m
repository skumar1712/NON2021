function [D,F] = FDmat(N,h)

%     Five points stencil
%     Both D and F are N x N matrices

    D11_ = [-25 48 -36 16 -3 zeros(1,N-5)]; 
    D1_ = [-3 -10 18 -6 1 zeros(1,N-5)] ; 
    D2_ = [1 -8 0 8 -1 zeros(1,N-5)]; 

    D1 = sparse(1:N-6,2:N-5, 1,N-6,N) ;
    D2 = sparse(1:N-6,3:N-4,-8,N-6,N) ;
%   D3 = sparse(1:N-6,4:N-3,0,N-6,N) ; 
    D4 = sparse(1:N-6,5:N-2,8,N-6,N) ; 
    D5 = sparse(1:N-6,6:N-1,-1,N-6,N) ; 
    D = D1 + D2 + D4 + D5; 
    Dn_ = [zeros(1,N-5) 1 -8 0 8 -1]; 
    Dn1_ = [zeros(1,N-5) -1 6 -18 10 3]; 
    Dnn1_ = [zeros(1,N-5) 3 -16 36 -48 25]; 
    D = [D11_; D1_ ; D2_ ; D ; Dn_ ; Dn1_; Dnn1_] /(12 * h) ;

    F11_ = [45 -154 214 -156 61 -10 zeros(1,N-6)] ; 
    F1_ = [10 -15 -4 14 -6 1 zeros(1,N-6)] ; 
    F2_ = [-1 16 -30 16 -1 zeros(1,N-5)] ; 

    F1 = sparse(1:N-6,2:N-5, -1,N-6,N) ;
    F2 = sparse(1:N-6,3:N-4,16,N-6,N) ;
    F3 = sparse(1:N-6,4:N-3,-30,N-6,N) ; 
    F4 = sparse(1:N-6,5:N-2,16,N-6,N) ; 
    F5 = sparse(1:N-6,6:N-1,-1,N-6,N) ; 
    F = F1 + F2 + F3 + F4 + F5; 
    
    Fn_ = [zeros(1,N-5) -1 16 -30 16 -1]; 
    Fn1_ = [zeros(1,N-6) 1 -6 14 -4 -15 10]; 
    Fn11_ = [zeros(1,N-6) -10 61 -156 214 -154 45]; 


    F = [F11_; F1_ ; F2_ ; F ; Fn_ ; Fn1_; Fn11_] /(12 * h^2) ;

return; 
