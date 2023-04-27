clc             % Clear command window
clear all       % Clear all variables in the workspace

% Author: Michal Ptak

R = 100;        % Resistance [ohm]
L = 0.04;       % Inductance [H]
C = 1*10^-5;    % Capacitance [F]
k = 0.9;        % Coefficient k [-]
M = k * sqrt(2*L^2); % Mutual inductance [H]

w = [           % Angular frequencies [rad/s]
    2*pi*50;
    3*2*pi*50;
    7*2*pi*50;
];

ZL = 1j*w.*L;   % Impedance of the inductors [ohm]
Zco = (-1j)./(w.*C); % Impedance of the capacitor [ohm]
ZM = 1j*w.*M;   % Impedance of the mutual inductor [ohm]

% Symbolic values of the circuit components
G = 1./R;       % Conductance [S]
YL = 1./ZL;     % Admittance of the inductors [S]
YC = 1./Zco;    % Admittance of the capacitor [S]
YM = 1./ZM;     % Admittance of the mutual inductor [S]

% Parameters of the sources
E1 = 100*sqrt(2);   % RMS value of the voltage of source 1 [V]
E2 = 70*sqrt(2);    % RMS value of the voltage of source 2 [V]
E3 = 50*sqrt(2);    % RMS value of the voltage of source 3 [V]

% Phasor voltages of the sources
Ea = [E1*exp(1j*0);
    E2*exp(1j*0);
    E3*exp(1j*0)];

Eb = [E1*exp(1j*(-2/3*pi));
    E2*exp(1j*0);
    E3*exp(1j*(2/3*pi))];

Ec = [E1*exp(1j*(2/3*pi));
    E2*exp(1j*0);
    E3*exp(1j*(-2/3*pi))];



% For loop to iterate over three phases
for i = 1:3; 

% Calculate impedances for each phase
ZA(i) = R + ZM(i);
ZB(i) = R + ZL(i) - ZM(i) + Zco(i);
ZC(i) = R + 2*ZL(i) - ZM(i) + (Zco(i)^2)/(2*Zco(i));

% Calculate admittances for each phase
YA(i) = 1./(ZA(i));
YB(i) = 1./(ZB(i));
YC(i) = 1./(ZC(i));

% Calculate the total admittance and the excitation vector for each phase
Y(i) = [YA(i)+YB(i)+YC(i)];
J(i) = [YA(i)*Ea(i)+YB(i)*Eb(i)+YC(i)*Ec(i)];

% Calculate the potential for each phase
V(i) = (Y(i))\(J(i));

% Calculate the branch currents for each phase
IA(i) = (Ea(i)-V(i))*YA(i);
IB(i) = (Eb(i)-V(i))*YB(i);
IC(i) = (Ec(i)-V(i))*YC(i);
 
% Calculate the voltages and powers for the loads in each phase

% Phase A
URA(i) = R*IA(i);
UMA(i) = ZM(i) * IA(i);

% Phase B
URB(i) = R * IB(i);
ULB(i) = ZL(i) * IB(i);
UMB(i) = -ZM(i) * IB(i);
UCB(i) = Zco(i) * IB(i);

% Phase C
URC(i) = R * IC(i);
ULC(i) = 2*ZL(i) * IC(i);
UCC(i) = (Zco(i)^2)/(2*Zco(i)) * IC(i);
UMC(i) = -ZM(i) * IC(i);

% Calculate the powers for the sources
SEA(i) = Ea(i)* conj(IA(i));
SEB(i) = Eb(i)* conj(IB(i));
SEC(i) = Ec(i)* conj(IC(i));

% Calculate the powers for the loads

% Phase A
SRA(i) = URA(i) * conj(IA(i));
SMA(i) = UMA(i) * conj(IA(i));
SA(i) = SRA(i)+SMA(i);

% Phase B
SRB(i) = URB(i) * conj(IB(i));
SLB(i) = ULB(i) * conj(IB(i));
SMB(i) = UMB(i) * conj(IB(i));
SCB(i) = UCB(i) * conj(IB(i));
SB(i) = SRB(i)+SLB(i)+SMB(i)+SCB(i);


% Phase C
SRC(i) = URC(i) * conj(IC(i));
SLC(i) = ULC(i) * conj(IC(i));
SMC(i) = UMC(i) * conj(IC(i));
SCC(i) = UCC(i) * conj(IC(i));
SC(i) = SRC(i)+SLC(i)+SMC(i)+SCC(i);
 end


% 1st HARMONIC BALANCE
BALANCE_W1 = IA(1)+IB(1)+IC(1)

% 3rd HARMONIC BALANCE
BALANCE_W3 = IA(2)+IB(2)+IC(2)

% 7th HARMONIC BALANCE
BALANCE_W7 = IA(3)+IB(3)+IC(3)

% VOLTAGE BALANCE
% 1st HARMONIC BALANCE
BALANCE_O11 = Ea(1) - Eb(1) - URA(1) - UMA(1) + UCB(1) + ULB(1) + UMB(1) + URB(1)
BALANCE_O12 = Eb(1) - Ec(1) - URB(1) - ULB(1) - UCB(1) - UMB(1) + UMC(1) + UCC(1) + ULC(1) + URC(1)

% 3rd HARMONIC BALANCE
BALANCE_O31 = Ea(2) - Eb(2) - URA(2) - UMA(2) + UCB(2) + ULB(2) + UMB(2) + URB(2)
BALANCE_O32 = Eb(2) - Ec(2) - URB(2) - ULB(2) - UCB(2) - UMB(2) + UMC(2) + UCC(2) + ULC(2) + URC(2)

% 7th HARMONIC BALANCE
BALANCE_O71 = Ea(3) - Eb(3) - URA(3) - UMA(3) + UCB(3) + ULB(3) + UMB(3) + URB(3)
BALANCE_O72 = Eb(3) - Ec(3) - URB(3) - ULB(3) - UCB(3) - UMB(3) + UMC(3) + UCC(3) + ULC(3) + URC(3)

% POWER BALANCE
% 1st HARMONIC BALANCE
POWER_BALANCE_1 = SEA(1) + SEB(1) + SEC(1) - SA(1) - SB(1) - SC(1)

% 3rd HARMONIC BALANCE
POWER_BALANCE_3 = SEA(2) + SEB(2) + SEC(2) - SA(2) - SB(2) - SC(2)

% 7th HARMONIC BALANCE
POWER_BALANCE_7 = SEA(3) + SEB(3) + SEC(3) - SA(3) - SB(3) - SC(3)

% PHASE CURRENT PLOT
IW_1 = [IA(1); IB(1); IC(1)];
IW_3 = [IA(2); IB(2); IC(2)];
IW_7 = [IA(3); IB(3); IC(3)];

% Plotting the current in the specified nodes for different harmonics using compass plot
figure (1);
subplot(1,3,1), compass(IW_1), title('Current in Node for 1st Harmonic')
subplot(1,3,2), compass(IW_3), title('Current in Node for 3rd Harmonic')
subplot(1,3,3), compass(IW_7), title('Current in Node for 7th Harmonic')

% Vector plot of voltages across the specified meshes for different harmonics
UO1_1 = [Ea(1); - Eb(1); - URA(1); - UMA(1); URB(1); ULB(1); UMB(1); UCB(1)];
UO1_2 = [Eb(1); - Ec(1); - URB(1); - ULB(1); - UMB(1); - UCB(1); URC(1); ULC(1); UMC(1); UCC(1)];

UO3_1 = [Ea(2); - Eb(2); - URA(2); - UMA(2); URB(2); ULB(2); UMB(2); UCB(2)];
UO3_2 = [Eb(2); - Ec(2); - URB(2); - ULB(2); - UMB(2); - UCB(2); URC(2); ULC(2); UMC(2); UCC(2)];

UO7_1 = [Ea(3); - Eb(3); - URA(3); - UMA(3); URB(3); ULB(3); UMB(3); UCB(3)];
UO7_2 = [Eb(3); - Ec(3); - URB(3); - ULB(3); - UMB(3); - UCB(3); URC(3); ULC(3); UMC(3); UCC(3)];

figure (2);
subplot(2,3,1), compass(UO1_1), title('Voltage Across Mesh 1 for 1st Harmonic')
subplot(2,3,4), compass(UO1_2), title('Voltage Across Mesh 2 for 1st Harmonic')

subplot(2,3,2), compass(UO3_1), title('Voltage Across Mesh 1 for 3rd Harmonic')
subplot(2,3,5), compass(UO3_2), title('Voltage Across Mesh 2 for 3rd Harmonic')

subplot(2,3,3), compass(UO7_1), title('Voltage Across Mesh 1 for 7th Harmonic')
subplot(2,3,6), compass(UO7_2), title('Voltage Across Mesh 2 for 7th Harmonic')