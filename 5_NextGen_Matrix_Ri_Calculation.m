%% 0. Next-Generation Matrix Construction

% Define symbolic variables for F
syms e S B1 B2 B3 B4 

% Matrix F definition
F = sym(zeros(6,6)); 
F(1,:) = [e*B3*S, e*B3*S, e*B4*S,e*B4*S, e*B1*S, e*B2*S];
F(2,:) = [(1-e)*B3*S, (1-e)*B3*S, (1-e)*B4*S, (1-e)*B4*S, (1-e)*B1*S, (1-e)*B2*S];

% Display the matrix F
disp(F);

%% 1. Matrix V Definition
syms v p_3 s_1 h_1 s_2 h_2 p_4 f_1

row1 = [p_3 + v, 0, 0, 0, 0, 0];
row2 = [0, p_3, 0, 0, 0, -(1 - s_1)*h_1];
row3 = [0, 0, p_3 + s_2*h_2 + v, 0, 0, 0];
row4 = [0, 0, 0, p_3 + s_2*h_2, 0, 0];
row5 = [-v, 0, 0, 0, p_4 + f_1, 0];
row6 = [0, 0, -v, 0, -f_1, (1 - s_1)*h_1 + s_1*h_1 + p_4];
V = [row1; row2; row3; row4; row5; row6];

V_invers = inv(V);
disp(V_invers);

%% 2. Next-Generation Matrix (Q)
Q = F * V_invers;

%eigenvalues
R_num = eig(Q);
disp(R_num);

%% 3. Substitute Parameters for Hospital Setting
syms B1 B2 B3 B4 S d g_1 g_2 x a l p_1 p_2 z
B1_new = d;
B2_new = d * (1 - g_1);
B3_new = d * x;
B4_new = d * x * (1 - g_2);
S_new = (l*z * (p_1+a) + a*l*(1-z)) / (p_2*(p_1+a));

% Replacement in the expression R_num function
R_final = subs(R_num, {B1, B2, B3, B4, S}, {B1_new, B2_new, B3_new, B4_new, S_new});
disp(R_final);

%% 4. Simplify Final Expression
R_final_simpl=simplify(R_final);
disp(R_final_simpl);
