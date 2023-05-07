%% Example for double compressor 
% 
clc
clear all;
%% Importing coolprop 
[v,e] = pyversion; 
system([e,' -m pip install --user -U CoolProp==6.3.0']);
import py.CoolProp.CoolProp.*
%% Degrees of Freedom 
Q_e=50000; %Required amount of refrigeration (W)
T_E=-20+273.15; %Evaporation temperature 
T_C=40+273.15; %Condensation temperature 
T_SH=-15+273.15; %Superheated vapor temperature 
T_SC=35+273.15; %Subcooled liquid temperature 
n_isen=0.85; %Isentropic compressor efficiency
n_me=0.90; %Mechanical efficiency 
beta=0.5; 
%% Refrigeration properties 
fluid='R134A';
P_E=PropsSI('P','T',T_E,'Q',0,fluid); %Evaporation pressure 
P_C=PropsSI('P','T',T_C,'Q',0,fluid); %Condensation pressure
P_I= (1-beta)*P_E+beta*P_C; %Intermediate pressure 
%% Defining the points 
H1=PropsSI('H','P',P_E,'T',T_SH,fluid); %First compressor inlet specific enthalpy
S1=PropsSI('S','P',P_E,'T',T_SH,fluid); %First compressor inlet specific entropy 
S2S=S1; %Ideal specific entropy at first compressor outlet 
H2S=PropsSI('H','P',P_I,'S',S2S,fluid); %Ideal specific entropy at first compressor outlet
H2=H1+(H2S-H1)/n_isen; %actual specific enthalpy at conderser inlet
H3=PropsSI('H','P',P_I,'Q',1,fluid); %intercooler gas specific enthalpy
H5=PropsSI('H','P',P_C,'T',T_SC,fluid); %Saturated vapor enthalpy at condensation temperature 
H6=H5; %Isoenthalpic expansion of first expansion valve 
H8=PropsSI('H','P',P_I,'Q',0,fluid); %intercooler liquid specific enthalpy
H9=H8; %Isoenthalpic expansion of second expansion valve
x6=(H6-H8)/(H3-H8); %vapor quality after first expansion valve 
H7= (1-x6)*H2+x6*H3; %Enthalpy after vapor mixture at intermidiate pressure 
S7=PropsSI('S','P',P_I,'H',H7,fluid); %Entropy at point 7
S4S=S7; %Ideal entropy at condensation inlet 
H4S=PropsSI('H','P',P_C,'S',S4S,fluid); %Ideal enthalpy at condensation inlet
H4=H7+(H4S-H7)/n_isen; %Actual enthalpy at condensation inlet
%% Saturation points at evaporation and condensation temperature
H10=PropsSI('H','P',P_E,'Q',1,fluid); %Specific enthalpy of saturated vapor at evaporation pressure
H11=PropsSI('H','P',P_C,'Q',1,fluid); %Specific enthalpy of saturated liquid at condensation pressure
H12=PropsSI('H','P',P_C,'Q',0,fluid); %Specific enthalpy of saturated liquid at condensation pressure
H13=PropsSI('H','P',P_E,'Q',0,fluid); %Specific enthalpy of saturated liquid at evaporation pressure
%% Required refrigerant mass flow rate 
m_e= Q_e/(H1-H9); %Required mass flow rate of evaporator 
m_c=m_e/(1-x6); %Condender mass flowrate 
m_i=m_c-m_e; %Saturated vapor mass flowrate at intermediate pressure 
%% Plotting Values 
P1=P_E; %Point 1
P2=P_I; %Point 2
P2S=P_I; 
P3=P_I; %Point 3
P4=P_C; %Point 4
P4S=P_C;
P5=P_C; %Point 5
P6=P_I; %Point 6
P7=P_I; %Point 7
P8=P_I; %Point 8
P9=P_E; %Point 9
P10=P_E; 
P11=P_C;
P12=P_C;
x9=(H9-H13)/(H10-H13); %Vapor quality after second expansion valve 
H_plot1=[H1 H2 H2S H7 H3 H6 H8 H9 H10 H1]; 
P_plot1=[P1 P2 P2S P7 P3 P6 P8 P9 P10 P1];
H_plot2=[H7 H4 H4S H11 H12 H5 H6 H3 H7]; 
P_plot2=[P7 P4 P4S P11 P12 P5 P6 P3 P7];
%% Creating function for S and T determination 
Svalues=@(x,y) PropsSI('S', 'P', x, 'H', y, fluid);
Tvalues=@(x,y) PropsSI('T', 'P', x, 'H', y, fluid);
 
S_plot1=[ Svalues(P1,H1) Svalues(P2,H2) Svalues(P2S,H2S) Svalues(P7,H7) Svalues(P3,H3) Svalues(P6,H6) Svalues(P8,H8) Svalues(P9,H9) Svalues(P10,H10) Svalues(P1,H1)]; 
T_plot1=[ Tvalues(P1,H1) Tvalues(P2,H2) Tvalues(P2S,H2S) Tvalues(P7,H7) Tvalues(P3,H3) Tvalues(P6,H6) Tvalues(P8,H8) Tvalues(P9,H9) Tvalues(P10,H10) Tvalues(P1,H1)]; 
S_plot2=[ Svalues(P7,H7) Svalues(P4,H4) Svalues(P4S,H4S) Svalues(P11,H11) Svalues(P12,H12) Svalues(P5,H5) Svalues(P6,H6) Svalues(P3,H3) Svalues(P7,H7)]; 
T_plot2=[ Tvalues(P7,H7) Tvalues(P4,H4) Tvalues(P4S,H4S) Tvalues(P11,H11) Tvalues(P12,H12) Tvalues(P5,H5) Tvalues(P6,H6) Tvalues(P3,H3) Tvalues(P7,H7)];  
%% Intial values 
P_Cr= 4000000; 
T_Cr=100+273.15; 
P_TL=[];
P_TV=[];
H_TL=[];
H_TV=[];
T_TL=[];
T_TV=[];
S_TL=[];
S_TV=[];
%%  T-s diagram
for T_loop = 273.15:1:T_Cr
    S_L = PropsSI('S', 'T', T_loop, 'Q', 0, fluid);
    S_TL=[S_TL S_L];
    T_TL=[T_TL T_loop];
    S_V = PropsSI('S', 'T', T_loop, 'Q', 1, fluid);
    S_TV=[S_V S_TV];
    T_TV=[T_loop T_TV];
end 
%%  P-h diagram
for P_loop = 4000:10000:P_Cr
    H_L = PropsSI('H', 'P', P_loop, 'Q', 0, fluid);
    P_L = PropsSI('P', 'P', P_loop, 'Q', 0, fluid);
    H_TL=[H_TL H_L];
    P_TL=[P_TL P_L];
    H_V = PropsSI('H', 'P', P_loop, 'Q', 1, fluid);
    P_V = PropsSI('P', 'P', P_loop, 'Q', 1, fluid);
    H_TV=[H_V H_TV];
    P_TV=[P_V P_TV];
end 
%% Plots 
figure 
% Plotting the P-h diagram 
plot(H_TL,P_TL) 
xlabel('Specific enthalpy (J/kg)')
ylabel('Pressure (Pa)')
title('P-h diagram of R134a')
hold on 
plot(H_TV,P_TV)
%Plotting the refrigeration cycle
labels1 = {'1','2','','7','3','6','8','9','',''};
plot (H_plot1, P_plot1,'-')
text(H_plot1,P_plot1,labels1)
labels2={'7', '4', '', '', '', '5', '6', '3', ''};
plot (H_plot2,P_plot2,'-')
text(H_plot2,P_plot2,labels2)
hold off 

figure 
% Plotting the T-s Diagram 
plot(S_TL,T_TL) 
xlabel('Specific entropy')
ylabel('Temperature (K)')
title('T-s diagram of R134a')
hold on 
plot(S_TV,T_TV)
% Plotting the refrigeration cycle 
plot (S_plot1, T_plot1,'-')
text(S_plot1,T_plot1,labels1)
plot (S_plot2, T_plot2,'-')
text(S_plot2,T_plot2,labels2)
hold off 
%% Printing Important values
fprintf('The properties of each points are given below:\n')
fprintf('%15s\t%8s\t%8s\t%8s\t%8s\n','Point','P','H','T','S')
fprintf('%15s\t%8s\t%8s\t%8s\t%8s\n','','(Pa)','(J/kg)','(K)','(J/K)')

fprintf('%15s\t%8.3f\t%8.3f\t%8.3e\t%8.3f\n','Point1',P1,H1,Tvalues(P1,H1),Svalues(P1,H1))
fprintf('%15s\t%8.3f\t%8.3f\t%8.3e\t%8.3f\n','Point2',P2,H2,Tvalues(P2,H2),Svalues(P2,H2))
fprintf('%15s\t%8.3f\t%8.3f\t%8.3e\t%8.3f\n','Point3',P3,H3,Tvalues(P3,H3),Svalues(P3,H3))
fprintf('%15s\t%8.3f\t%8.3f\t%8.3e\t%8.3f\n','Point4',P4,H4,Tvalues(P4,H4),Svalues(P4,H4))
fprintf('%15s\t%8.3f\t%8.3f\t%8.3e\t%8.3f\n','Point5',P5,H5,Tvalues(P5,H5),Svalues(P5,H5))
fprintf('%15s\t%8.3f\t%8.3f\t%8.3e\t%8.3f\n','Point6',P6,H6,Tvalues(P6,H6),Svalues(P6,H6))
fprintf('%15s\t%8.3f\t%8.3f\t%8.3e\t%8.3f\n','Point7',P7,H7,Tvalues(P7,H7),Svalues(P7,H7))
fprintf('%15s\t%8.3f\t%8.3f\t%8.3e\t%8.3f\n','Point8',P8,H8,Tvalues(P8,H8),Svalues(P8,H8))
fprintf('%15s\t%8.3f\t%8.3f\t%8.3e\t%8.3f\n','Point9',P9,H9,Tvalues(P9,H9),Svalues(P9,H9))
fprintf ('\nThe required R134a evaporator mass flowrate is=%f kg/s\n', m_e);
%% Energy Balance
Q_eva= m_e*(H1-H9); %Evaporation rate 
Q_cond=(m_e/(1-x6))*(H4-H5); %Condensation rate 
W_comp=m_e*(H2-H1)+(m_e/(1-x6))*(H4-H7); %Compressor power 
Energy_balance=Q_eva-Q_cond+W_comp; 
if abs (Energy_balance)<1 %Setting a limit of error 
    fprintf('\nThe energy balance verifies the solution\n')
else 
    fprintf ('\nError!! The energy is not balanced, check the solution again.\n')
end 