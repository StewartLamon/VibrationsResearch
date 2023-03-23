%Timoshenko Version Stepped Beam
clear all , close all, clc

E=[1,1,1]; %Modulus of elasticity for each segment, Pa
L=[1]; %Length of beam, 10 in. (m)
p=[1, 1, 1]; %Density of segments (kg/m^3)
W = [1, 1, 1]; %beam width, 1 inch (m)
H = [1, 1, 1]; %beam height, quarter of an inch (m)
h = [0]; %void cross section size array (m)
l = 0.0127; %1/2 in void length (m)
A=[W(1)*H(1), W(2)*H(2), W(3)*H(3)]; %Area of each segment
I=[1/12*W(1)*H(1)^3, 1/12*W(2)*H(2)^3 1/12*W(3)*H(3)^3]; %Second moment area for each segment
e=[0, 0.5]%,0.09101582]; %Dimensionless spatial coordinate for each discontinuity
segments = length(e); %Number of segments of beam based on discontinuities

%Timoshenko specific variables
V = 1; %Poisson's Ratio
K = 1; % Shear correction factor
G = [1, 1, 1]; % Shear Modulus


%Setting first segment lambda as variable
syms lambda_1

lambda(1) = lambda_1;

%Evaluating all lambdas in terms of first segment lambda
for i=2:segments
    lambda(i)=lambda_1*((E(1)*I(1))/(p(1)*A(1)))*((p(i)*A(i))/(E(i)*I(i)));
end

%Defining other terms for Timoshenko (these are both squared):

%Dimensionless r squared
for i=1:segments
    r(i) = I(i)/(A(i)*L^2);
end

%Dimensionless s squared
for i = 1:segments
    s(i) = E(i)*I(i)/(K*G(i)*A(i)*L^2);
end

%Dimensionless d variables

d_one = lambda_1^4*(r(1) + s(1))/2;

for i = 1:segments
    d_two(i) = (lambda(i)^4)*((lambda(i)^4)*s(i)*r(i) - 1);
end

%Defining Beta variables

for i = 1:segments
    Beta_one(i) = (-d_one + (d_one^2 - d_two(i))^(1/2))^(1/2);
end

for i = 1:segments
    Beta_two(i) = (d_one + (d_one^2 - d_two(i))^(1/2))^(1/2);
end

%Defining M variables
for i = 1:segments
    m_one(i) = (Beta_one(i)^2 + (lambda(i)^4)*s(i))/Beta_one(i);
end

for i = 1:segments
    m_two(i) = (-Beta_two(i)^2 + (lambda(i)^4)*s(i))/Beta_two(i);
end

%Defining mu

for i = 1:segments-1
    mu(i) = A(i+1)/A(i);
end

%Coefficient for continuity between segments
for i=1:segments-1
    Xi(i)=E(i+1)*I(i+1)/(E(i)*I(i)); % assign value to the ith element
end

%Creating P matrix
for i=1:segments-1
    P(:,:,i)=[1,0,1,0;
        0, m_one(i+1), 0, -m_two(i+1);
        m_one(i+1)*Beta_one(i+1)*Xi(i), 0, m_two(i+1)*Beta_two(i+1)*Xi(i), 0;
        0, (m_one(i+1) - Beta_one(i+1))*mu(i), 0, -(m_two(i+1) - Beta_two(i+1))*mu(i)*Xi(i)];
end

%delta=difference between locations of segments i and i-1
for i = 2:segments
    delta(i-1) = e(i) - e(i-1);
end

%Creating Q matrix
for i=1:segments-1
    Q(:,:,i)=[cosh(Beta_one(i)*delta(i)),  sinh(Beta_one(i)*delta(i)), cos(Beta_two(i)*delta(i)), sin(Beta_two(i)*delta(i));
        m_one(i)*sinh(Beta_one(i)*delta(i)), m_one(i)*cosh(Beta_one(i)*delta(i)),  m_two(i)*sin(Beta_two(i)*delta(i)),  -m_two(i)*cos(Beta_two(i)*delta(i));
        m_one(i)*Beta_one(i)*cosh(Beta_one(i)*delta(i)),  m_one(i)*Beta_one(i)*sinh(Beta_one(i)*delta(i)),  m_two(i)*Beta_two(i)*cos(Beta_two(i)*delta(i)),  m_two(i)*Beta_two(i)*sin(Beta_two(i)*delta(i));
        (m_one(i) - Beta_one(i))*sinh(Beta_one(i)*delta(i)),   (m_one(i) - Beta_one(i))*cosh(Beta_one(i)*delta(i)),  (m_two(i) + Beta_two(i))*sin(Beta_one(i)*delta(i)),  -(m_two(i) + Beta_two(i))*cos(Beta_one(i)*delta(i))];
end

%Applying Transfer Matrix
for i=1:segments-1
    T(:,:,i)=inv(P(:,:,i))*Q(:,:,i);
end
fprintf('Finished inversion \n')
alpha=1-e(segments);
%Clamped boundary condition

Capital_Lambda=[m_one(segments)*Beta_one(segments)*cosh(Beta_one(segments)*alpha), m_one(segments)*Beta_one(segments)*sinh(Beta_one(segments)*alpha), m_two(segments)*Beta_two(segments)*cos(Beta_two(segments)*alpha), m_two(segments)*Beta_two(segments)*sin(Beta_two(segments)*alpha);
   (m_one(segments) - Beta_one(segments))*sinh(Beta_one(segments)*alpha), (m_one(segments) - Beta_one(segments))*cosh(Beta_one(segments)*alpha), (m_two(segments) + Beta_two(segments))*sin(Beta_two(segments)*alpha), -(m_two(segments) - Beta_two(segments))*cos(Beta_two(segments)*alpha)];

%Loop to define Capital Gamma
C_Gam = Capital_Lambda;
for n=1:segments-1
    C_Gam=C_Gam*T(:,:,segments-n);
end
%
fprintf('Finished Cap_Gam \n')
% m_one_final = 1;
% m_two_final = 1;

%Creating factored m matrices?
% for i = 1:segments
%     m_one_final = m_one_final * m_one(i);
%     m_two_final = m_two_final * m_two(i);
% end
%System of equations representing disp, shear, etc relationships between
%all segments

system=[C_Gam;
    1,0,1,0;
    0, m_one(1), 0 , -m_two(1)];
fprintf('Made System \n')
%Simplifying and solving equation

%Not loading here
equation=det(system)
fprintf('Made Equation \n')
%
% Initial guesses (from book)
N = [1.8751, 4.69409, 7.854757, 10.995540, 14.137168, 17.2787595];

%simplify(equation)
fun=matlabFunction(equation)

for i = 1:6
    %digits(1)
    lam(i) = fsolve(fun,N(i));
    fprintf('Finished frequency %4.5f \n',i')
end

Frequencies = (lam.^4*E(1)*I(1)./(p(1)*A(1)*L.^4)).^((1/2))'
%