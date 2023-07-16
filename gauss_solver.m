%
% "Program for solving systems of equations using the Gaussian elimination method
% 
close all;
clear all;
n=input('Enter the number of variables  : ')
Rw=zeros(n,n+1)
for i=1:n
    for k=1:n+1
    Rw(i,k) =input('Enter the next parameter of the equation : ')
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
A = zeros(1,n+1); 
for p=1:n
    rn=n-p+1;
    minus=0;
    for s=1:n           
        minus= minus +(A(s)*Rw(rn,s));
        
    end
    A(rn) =(Rw(rn,n+1)-minus)/Rw(rn,rn);
end
Rw
for d=1:n
X=['X_',num2str(d),' is equal: ',num2str(A(d))];
disp(X)
end


