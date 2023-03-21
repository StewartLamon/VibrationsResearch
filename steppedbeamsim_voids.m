clear all , close all, clc

segments=3; %number of segments

cv2m=0.0254;%converts in to meters

E=20*10^10; %modulus
L=10*cv2m; %Total Length m
p=7850; %density kg/m^3
e=[0,1/3-.25*cv2m/L,1/3+0.25*cv2m/L]; %location of the begining of section i in dimensionless coordinates

iterations=0:1/64*cv2m:4/16*cv2m;
for loop=1:length(iterations)

    %redefine the area and moment of inertia in terms of the base and height
    b=1*cv2m;
    h=.25*cv2m;

    for i=1:segments
        I(i)=b*h^3/12;
        A(i)=b*h;
    end

    %void moment and area adjustments
    voidlocation=2; %void is at section i
    bv=0.1*cv2m;%2/16*cv2m;
    hv=iterations(loop);%4/16*cv2m;

    for i=voidlocation
        I(i)=I(i)-bv*hv^3/12;
        A(i)=A(i)-bv*hv
    end

    %write all lambdas in terms of lambda_1
    syms lambda_1

    lambda(1)=lambda_1;

    for i=2:segments
        lambda(i)=lambda_1*(I(1)/A(1))*(A(i)/I(i));
    end

    %calculate all Xi
    for i=1:segments-1
        Xi(i)=I(i+1)/I(i);
    end

    %Calculate all P matricies
    for i=1:segments-1
        P(:,:,i)=[1,0,1,0;
            0, lambda(i+1), 0, lambda(i+1);
            (lambda(i+1)^2)*Xi(i), 0, -(lambda(i+1)^2)*Xi(i), 0;
            0, (lambda(i+1)^3)*Xi(i), 0, -(lambda(i+1)^3)*Xi(i)];
    end

    %delta=difference between locations of segments i and i-1
    %can also be thought of as the length of segment i
    for i=1:segments-1
        delta(i)=e(i+1)-e(i);
    end

    %Calculate all Q Matricies
    for i=1:segments-1
        Q(:,:,i)=[cosh(lambda(i)*delta(i)),  sinh(lambda(i)*delta(i)), cos(lambda(i)*delta(i)), sin(lambda(i)*delta(i));
            lambda(i)*sinh(lambda(i)*delta(i)),  lambda(i)*cosh(lambda(i)*delta(i)),  -lambda(i)*sin(lambda(i)*delta(i)),  lambda(i)*cos(lambda(i)*delta(i));
            (lambda(i)^2)*cosh(lambda(i)*delta(i)),  (lambda(i)^2)*sinh(lambda(i)*delta(i)),  -(lambda(i)^2)*cos(lambda(i)*delta(i)),  -(lambda(i)^2)*sin(lambda(i)*delta(i));
            (lambda(i)^3)*sinh(lambda(i)*delta(i)),  (lambda(i)^3)*cosh(lambda(i)*delta(i)),  (lambda(i)^3)*sin(lambda(i)*delta(i)),  -(lambda(i)^3)*cos(lambda(i)*delta(i))];
    end

    %Calculate All T Matricies
    %This relates the A,B,C,D constants of the ith section to the i+1 section
    for i=1:segments-1
        T(:,:,i)=inv(P(:,:,i))*Q(:,:,i);
    end

    alpha=1-e(segments);
    %Right Side free boundary conditions
    Capital_Lambda=[cosh(lambda(segments)*alpha),sinh(lambda(segments)*alpha),-cos(lambda(segments)*alpha),-sin(lambda(segments)*alpha);
        sinh(lambda(segments)*alpha),cosh(lambda(segments)*alpha),sin(lambda(segments)*alpha),-cos(lambda(segments)*alpha)];

    %loop to define Capital Gamma
    C_Gam=Capital_Lambda;
    for n=1:segments-1
        C_Gam=C_Gam*T(:,:,segments-n);
    end

    %combination of left and right side boundary conditions
    %left side clamped, right side free
    system=[C_Gam(1,1),C_Gam(1,2),C_Gam(1,3),C_Gam(1,4);
        C_Gam(2,1),C_Gam(2,2),C_Gam(2,3),C_Gam(2,4);
        1,0,1,0;
        0,1,0,1];

    %to solve the system for nontrivial zeros the determinant has to be zero
    equation=det(system)

    simplify(equation)

    %calculate all dimensionless eigenvalues up to ___
    n=1;
    lambda_s(1,loop)=0;
    error=0.01;
    for i=1:15%4.5:0.01:5
        check=vpasolve(equation,lambda_1,i);
        if(~(check < lambda_s(n,loop)+error) && (check > lambda_s(n,loop)-error))
            lambda_s(n+1,loop)=check
            n=n+1;
        end
    end

    %calculate all corresponding time and space frequencies omega
    for n=1:length(lambda_s(:,loop))
        omega(n,loop)=sqrt(((lambda_s(n,loop))^(4))*E*I(1)/(p*A(1)*L^4))
        beta_1(n,loop)=lambda_s(n,loop)*L;
    end

end

%need to find all of the constants A,B,C,D
