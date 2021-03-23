clear
close all

Batch_no = 1; % batch number =1
Num_of_Batches = 1;
%% par

intial_conds = 0.5;

x0.mux = 0.41;
x0.mup = 0.041;
x0.S = 1;
x0.DO2 = 15;
x0.X = intial_conds;
x0.P = 0;
x0.V =  5.800e+04;
x0.Wt = 6.2e+04;
x0.CO2outgas = 0.038;
x0.O2 = 0.20;
x0.pH = 6.5;
x0.T = 297;
x0.a0 = intial_conds*(1/3);
x0.a1 = intial_conds*(2/3);
x0.a3 = 0;
x0.a4 = 0;
x0.Culture_age = 0;
x0.PAA = 1400;
x0.NH3 = 1700;
alpha_kla = 85; 
PAA_c = 530000;
N_conc_paa = 2*75000; 

h = 0.2;
par = Parameter_list(x0,alpha_kla,N_conc_paa,PAA_c);
T =230;
Batch_time = 0:h:T;

%%  Set up simulation flags
Ctrl_flags.SBC= 0;      % ignore   

Batch_run_flags.Control_strategy = ones(1,Num_of_Batches);
Batch_run_flags.Batch_length = 0*ones(1,Num_of_Batches);
Batch_run_flags.Batch_fault_order_reference= 0*ones(Num_of_Batches,1);
Batch_run_flags.Raman_spec = 0*ones(1,Num_of_Batches);

Ctrl_flags.PRBS= Batch_run_flags.Control_strategy(Batch_no);    % 0 - Recipe driven (i.e Sequential batch control (SBC)) 
                        % 1- Operator conntroller batches
Ctrl_flags.Fixed_Batch_length = Batch_run_flags.Batch_length(Batch_no); % 0 - Fixed Batch length
                                   % 1 - Uneven batch length
Ctrl_flags.IC = 0;      % Initial Conditions flag (IC) = 0 -  Randomly calculated initial conditions must control using SBC
Ctrl_flags.Inhib = 2;   % Inhibition  flag (Inhib)     = 0 - No inhibition
                        % Inhibition  flag (Inhib)     = 1 - Inhibition DO2, T, pH
                        % Inhibition  flag (Inhib)     = 2 - Inhibition of DO2,T,pH,CO_2_L,PAA and N
Ctrl_flags.Dis = 1;     % Disturbance flag  (Dis)      = 0 - No process disturbances
                        % Disturbance flag  (Dis)      = 1 - In batch fluctuations on mu_P, mu_x, c_s, c_oil,abc,PAA_c,Tcin and O2_in            
Ctrl_flags.Faults= Batch_run_flags.Batch_fault_order_reference(Batch_no);   % Fault flag  (Faults)         = 0 - No Faults
                        % Fault flag  (Faults)         = 1 - Aeration rate fault
                        % Fault flag  (Faults)         = 2 - Vessel back pressure  fault
                        % Fault flag  (Faults)         = 3 - Substrate feed rate fault
                        % Fault flag  (Faults)         = 4 - Base flowrate fault
                        % Fault flag  (Faults)         = 5 - Coolant flowrate fault
                        % Fault flag  (Faults)         = 6 - All of the above faults
                        % Fault flag  (Faults)         = 7 - Temperature sensor error
                        % Fault flag  (Faults)         = 8 - pH sensor error
Ctrl_flags.Vis= 0;      % Viscosity flag (Vis)          = 0 - uses simulated viscosity
Ctrl_flags.Raman_spec = Batch_run_flags.Raman_spec(Batch_no);  %  Raman_spec =  0 - No spectral data recorded
                            %  Raman_spec =  1 - Spectral data recorded
                            %  Raman_spec =  2 - Spectral data used to  predict and control PAA
Ctrl_flags.Batch_Num =  Batch_no;                     

                                               
% Off-line measurement sampling rate and analysis delay
Ctrl_flags.Off_line_m =  12;     % Off-line measurement sampling rate (hours)
Ctrl_flags.Off_line_delay =  4;  % Off-line measurement analysis time delay (hours)
% Plot graphs 
Ctrl_flags.plots =  1;           % 0 - No plots 1 - plots

%% Temperature and pH Set points for Batch
Ctrl_flags.T_sp = 298; % Temperature Set-point (K)
Ctrl_flags.pH_sp = 6.5; % pH Set-point (-)
    
%% Creates process disturbances on growth rates as well as process inputs
% using a low pass filter
b1 = 1 - 0.995;
a1 = [1 -0.995];
% Penicillin specific growth rate disturbance: with SD of +/- 0.0009 [hr^{-1}]
v = zeros(T/h+1, 1);
distMuP = filter(b1,a1,0.03*v);
Xinterp.distMuP = createChannel('Penicillin specific growth rate disturbance','g/Lh','h',Batch_time,distMuP);
% Biomass specific growth rate disturbance: with SD  +/- 0.011 [hr^{-1}]
% v = randn(T/h+1, 1);
distMuX = filter(b1,a1,0.25*v);
Xinterp.distMuX = createChannel('Biomass specific  growth rate disturbance','hr^{-1}','h',Batch_time,distMuX);
% Substrate inlet concentration disturbance: +/- 15 [g L^{-1}]
% v = randn(T/h+1, 1);
distcs = filter(b1,a1,5*300*v);
Xinterp.distcs = createChannel('Substrate concentration disturbance ',' g L^{-1}','h',Batch_time,distcs);
% Oil inlet concentration disturbance: +/- 15 [g L^{-1}]
% v = randn(T/h+1, 1);
distcoil = filter(b1,a1,300*v);
Xinterp.distcoil = createChannel('Substrate concentration disturbance ',' g L^{-1}','h',Batch_time,distcoil);
% Acid/Base molar inlet concentration disturbance: +/- 0.004 [mol L^{-1}]
% v = randn(T/h+1, 1);
% distabc = filter(b1,a1,0.2*v);
distabc = filter(b1,a1,0.2*v);
Xinterp.distabc = createChannel('Acid/Base concentration disturbance ',' g L^{-1}','h',Batch_time,distabc);
% v = randn(T/h+1, 1);
distPAA = filter(b1,a1,300000*v);
Xinterp.distPAA = createChannel('Phenylacetic acid concentration disturbance ',' g L^{-1}','h',Batch_time,distPAA);
% PAA inlet concentration disturbance: +/- 120  [g  L^{-1}]
% v = randn(T/h+1, 1);
distPAA = filter(b1,a1,300000*v);
Xinterp.distPAA = createChannel('Phenylacetic acid concentration disturbance ',' g L^{-1}','h',Batch_time,distPAA);
% Coolant temperature inlet concentration disturbance: +/- 3 [K]
% v = randn(T/h+1, 1);
% distTcin = filter(b1,a1,100*v);
distTcin = filter(b1,a1,100*v);
Xinterp.distTcin = createChannel('Coolant inlet temperature disturbance ','K','h',Batch_time,distTcin);
% Oxygen inlet concentration: +/- 0.009 [%]
% v = randn(T/h+1, 1);
distO_2in = filter(b1,a1,0.02*v);
Xinterp.distO_2in = createChannel('Oxygen inlet concentration','%','h',Batch_time,distO_2in);



%% run
[Xref] = indpensim(@fctrl_indpensim, Xinterp, x0, h, T,2,par,Ctrl_flags);

%% plot
All_variables_names =fieldnames(Xref);
Var_all = size(All_variables_names, 1);
Batch_start = 1;
P_finalvalve = Xref.P.y(end);