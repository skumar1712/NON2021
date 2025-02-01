% -------------------------------------------------------
% The following function constructs the initial data i.e. the hyperbolic
% planar polygon
% Arguments: 
% M: Total number of sides
% N: The discretized nodes
% l = L /(M-1)
% -------------------------------------------------------

function [X,T] = GetPol(N,M,l,b)

% ---------------------------
% The tangent vectors
% ---------------------------
    a = sqrt(1+b^2) ;
    k = 0 : M-1 ; k=k.' ;
    T_vert = [a*cosh((k+(1-M)/2)*l) a*sinh((k+(1-M)/2)*l) b*ones(M,1)];
    
    for k = 0 : M-1
        T(k * N/M + (1:N/M),:) = ones(N/M,1) * T_vert(k+1,:) ; 
    end 

% -----------------------------------------------------
% drawing the polygon i.e. by integrating X_s = T 
% -----------------------------------------------------  
  
    L = M*l ; h = L / N ; 

    % ----- 
    % Computing X from T
    % ----- 
%     Xc1 = h*cumsum(T); Xc1 = [0 0 0 ; Xc1] ; 
%     X1mean = h*0.5* (Xc1(1,1) + 2*sum(Xc1(2:end-1,1)) + Xc1(end,1))/L ;
%     X2mean = h*0.5* (Xc1(1,2) + 2*sum(Xc1(2:end-1,2)) + Xc1(end,2))/L ;
%     X3mean = h*0.5* (Xc1(1,3) + 2*sum(Xc1(2:end-1,3)) + Xc1(end,3))/L ;
%     Xc1 = Xc1- [X1mean X2mean X3mean] ; 

% ----- 
    % Computing X explicitly
    % -----   
%     k = -M/2 : -M/2+1;  
%     X12 = (l/2) * a  * [sinh(k*l) ; cosh(k*l)] / (sinh(l/2)) ;
%     Xx1 = linspace(X12(1,1),X12(1,2),N/M+1); 
%     Xx2 = linspace(X12(2,1),X12(2,2),N/M+1); 
%     Xx = [Xx1 ; Xx2]; 
%     R = [cosh(l) sinh(l); sinh(l) cosh(l)] ; 
%     for j = 0 : M-2
%         XXaux = R^j*Xx ;
%         XX(:,j*N/M+(1:N/M)) = XXaux(:,1:N/M);
%     end
%     j = M-1; 
%     XXaux = R^j*Xx ;
%     XX(:,j*N/M+(1:N/M+1)) = XXaux;
    
    k = -M/2 : -M/2+1;  
    X12 = (l/2) * a  * [sinh(k*l) ; cosh(k*l)] / (sinh(l/2)) ;
    Xx1 = linspace(X12(1,1),X12(1,2),N/M+1); 
    Xx2 = linspace(X12(2,1),X12(2,2),N/M+1); 
    Xx = [Xx1 ; Xx2]; 
    R = [cosh(l) sinh(l); sinh(l) cosh(l)] ; 
    for j = 0 : (M-2)/2-1
        XXaux = R^j*Xx ;
        XX(:,j*N/M+(1:N/M)) = XXaux(:,1:N/M);
    end
    j = (M-2)/2; 
    XXaux = R^j*Xx ;
    XX(:,j*N/M+(1:N/M+1)) = XXaux;
    XX = XX.' ; 
    
    XXX = [XX(:,1) XX(:,2); -flipud(XX(1:end-1,1)) flipud(XX(1:end-1,2))];
    
    s = linspace(-L/2,L/2,N+1); 
    X = [XXX b*s.' ] ; 
    % -------------------------------------------------------
%    Making T a N+1 point vector by taking the average of points

    TT = [T(1, :); .5 * (T(1:end-1, :) + T(2:end, :)); T(end, :)] ; 

    T = [TT(:, 1)  TT(:, 2) TT(:,3)] ./ sqrt(TT(:, 1).^2 - TT(:, 2).^2- TT(:, 3).^2) ; 
    
end


