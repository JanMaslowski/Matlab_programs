%
% Program for function interpolation using Czebyszow polynomials.
%
close all;
clear all;
% Enter 3 or more points
n=input('Enter the number of points : ')
Rw=zeros(n,n+1)
for i=1:n
    vector=input('point coordinates E.g. "[1,1]"   : ')
    for k=3:n
     
     Rw(i,1)=1;
     Rw(i,2)=vector(1);
     Rw(i,n+1)=vector(2);
     Rw(i,k)=(2*vector(1)*Rw(i,k-1)-Rw(i,k-2))
    end
end
sw=Rw;
swsquare=sw(:, 1:n);
if det(swsquare)==0
    disp('The determinant of the matrix is equal to 0, the system has no solution.')
    return
end
for z=1:n-1
    for k=z+1:n
        count=1;
        for i=1:n
            j=n+1-i;
            Rw(k,j)=Rw(k,j)-(Rw(k,z)./Rw(z,z)).*Rw(z,j);
                if count==1
                Rw(k,n+1)=Rw(k,n+1)-(Rw(k,z)./Rw(z,z)).*Rw(z,n+1);
                end
                count=count+1;
        end
        
    end
end
A = zeros(1,n); 
for p=1:n
    rn=n-p+1;
    minus=0;
    for s=1:n           
        minus= minus +(A(s)*Rw(rn,s));
        end
    A(rn) =(Rw(rn,n+1)-minus)/Rw(rn,rn);
end
A;
disp('Parameters of the interpolated function are: .')
for d=1:n
X=['X_',num2str(d),' is equal: ',num2str(A(d))];
disp(X)
end
hold on
x=linspace(0.9,5);
syms sx;
syms W
for k=3:n
     W(1)=1;
     W(2)=sx;
     W(k)=((2.*sx.*W(k-1))-W(k-2));
end
F=A*W';
F_numeric=subs(F,sx,x);
plot(x,F_numeric);
for r = 1:n
  scatter(sw(r,2),sw(r,n+1),'filled',MarkerFaceColor='r');  
end
title('Interpolation with Chebyshev polynomials.')