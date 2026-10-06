function [iter,SOL,incre,I,tcl,tcc,incre1,ACOC,dif_iter] = Newton_matrcial(A,X0,maxiter,tol)
%[iter,SOL,incre,I,tcl,tcc,incre1,ACOC,dif_iter] = Newton_matrcial(A,X0,maxiter,tol)
n = size(A); n = n(1); 
incre1 = norm(eye(n)-A*X0);
I = []; %vector de incrementos. 
I = [I,incre1];
iter = 0; 
incre = incre1; 
dif_iter = []; %diferencia entre dos iterados consecutivos
%Utilizmaos como incremento I-A*X0 pues la sol es la inversa de A
if (incre<1)
    while(incre>tol && iter<maxiter)
        X = X0*(2*eye(n)-A*X0); 
        incre = norm(eye(n)-A*X0); 
        Dif = norm(X-X0); 
        X0 = X; 
        iter = iter+1; 
        I = [I,incre];
        dif_iter = [dif_iter,Dif];
    end 
    
    if (iter >= maxiter)
        disp('Se ha superado el número máximo de iteraciones') 
        SOL=[]; 
    else 
        SOL = X;  

    end 
    figure; 
    ACOC = log(I(3:end)./I(2:end-1))./log(I(2:end-1)./I(1:end-2)); 
    tcl=I(2:end)./I(1:end-1);subplot(3,1,1);plot(tcl);legend('Tasa lineal')
    tcc=I(2:end)./(I(1:end-1).^2);subplot(3,1,2);plot(tcc);legend('Tasa cuadrática')
    subplot(3,1,3);
    plot(ACOC);
    legend('ACOC')
else 
    disp('No se cumple la condición para la norma del error')
    SOL=[];
end 
