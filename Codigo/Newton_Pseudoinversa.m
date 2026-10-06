function [SOL, iter, I, incre, tcl,tcc,ACOC,dif_iter] = Newton_pseudoinversa(A,X0,maxiter,tol)
%[SOL, iter, I, incre, tcl,tcc,ACOC,dif_iter] = Newton_pseudoinversa(A,X0, maxiter,tol)
t = size(X0);
n = t(2);

incre1(1) = norm(A*X0*A-A);
incre1(2) = norm((A*X0)'-A*X0); 
incre1(3) = norm(X0*A*X0-X0);
incre1(4) = norm((X0*A)'-X0*A); 
incre = sum(incre1);
iter = 1; 
I = []; I=[I,incre]; 
dif_iter = []; 
while(iter<maxiter && incre > tol)
    XK = X0*(2*eye(n)-A*X0);
    dif = norm(XK-X0); 
    incre1(1) = norm(A*XK*A-A);
    incre1(2) = norm((A*XK)'-A*XK); 
    incre1(3) = norm(XK*A*XK-XK); 
    incre1(4) = norm((XK*A)'-XK*A); 
    incre = sum(incre1);
    iter = iter+1;
    X0 = XK; 
    I=[I,incre];
    dif_iter = [dif_iter,dif]; 
end  
if iter >= maxiter 
    disp('Se ha superado el numero maximo de iteraciones')
    SOL = []; 
else 
    disp('Se ha encontrado la solución')
    SOL = XK; 
end 
    figure; 
    ACOC = log(I(3:end)./I(2:end-1))./log(I(2:end-1)./I(1:end-2)); 
    tcl=I(2:end)./I(1:end-1);
    tcc=I(2:end)./(I(1:end-1).^2);
end 