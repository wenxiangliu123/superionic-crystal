clear all
% Setting for correlations
% Ns = 1;         % sample interval in lammps
Nc_all = 500;       % number of data for correlation
Nm_all = 500;
Ncm_all = 500;
total_time_step = 50000;

%% settings for MD
time_step = 0.001;                                            % timestep in MD
kB = 8.617343e-5;                                             % Boltznman, unit is eV/K
ev2J = 1.6022e-19;
ps2s = 1e-12;
a2m = 1e-10;
g2kg = 1e-3;
NA = 6.022e23;

%% read data and begin to calculate

finish = 'Begin to read log.lammps';
%read file
maindir = 'D:\program2-ionic_crystal\ML-potential\Ag2Te\LJ-test\entralpy\0.5-0.05\mass\107-30\800';
subdir =  dir( maindir );       % the name of sub file  
size0 = size(subdir);
length = size0(1);
c1 = 612;                       % lines reqiured to skip
T1 = 10;                       % total number of lines are read
Nlog1 = 6;                     % colunmns to be read

V = [];                         % cell length
cell_x = [];cell_y = [];cell_z = [];
Temperature = [];               % temperature against time
NN = length-2;                  % numbers of simulations

hac_all = [];hac_cond=[]; hac_ct = [];hac_pt = [];
mac_all=[];hmac_all=[];
mace_all = []; hmace_cond = [];hmace_cond_2 = [];hmac_cond = [];
hac_cpt = [];hac_pct = [];

%hac_fft = [];mac_fft=[];hmac_fft=[];
%H_Ag = []; H_Te=[];
    
for ii= 3:length

    % read lammps.log
    Formatstring1 = repmat('%.16f', 1, Nlog1); % format for lammps.log
    file1 = fullfile( maindir, subdir(ii).name, 'log.lammps' );
    mFID = fopen(file1); 
    y = textscan( mFID, Formatstring1, T1, 'Headerlines', c1);
    data(1: T1, :) = cell2mat(y);

    Temperature = [Temperature, data(:, 2)];
    V = [V, data(1, Nlog1-2)*data(1, Nlog1-1)*data(1, Nlog1)];
    %cell_length = [data(1, 4)+data(1, 5)+data(1, 6)]/3;
    cell_x = [cell_x, data(1, Nlog1-2)];
    cell_y = [cell_y, data(1, Nlog1-1)];
    cell_z = [cell_z, data(1, Nlog1)];
    clear data;
    fclose('all');
    finish = 'Begin to read heat.txt';

    % read heat current file
    file_all_1 = fullfile( maindir, subdir(ii).name, 'heat_momentum_flux.txt' );
    %[Time_all, Jcx_all, Jpx_all, Jcy_all, Jpy_all, Jcz_all, Jpz_all,Jx_mass,Jy_mass,Jz_mass]...
    %= textread(file_all, '%n%n%n%n%n%n%n%n%n%n', 'headerlines', 2);

    [Time_all, Jx_all, Jy_all, Jz_all, Jcx_all, Jcy_all, Jcz_all, Jx_mass,Jy_mass,Jz_mass]...
    = textread(file_all_1, '%n%n%n%n%n%n%n%n%n%n', 'headerlines', 2); 

    file_all_2 = fullfile( maindir, subdir(ii).name, 'enthalpy.txt' );
    [~, h_Ag, h_Te]...
    = textread(file_all_2, '%n%n%n', 'headerlines', 2);
    
    %H_Ag = [H_Ag;h_Ag];H_Te = [H_Te;h_Te];
    
    %data_mass_flux = load([maindir, '\',subdir(ii).name,'\','Ag_mass_flux.mat']);
    
    %J_mass = data_mass_flux.mass_flux;
    %Jx_mass = J_mass(:,1,:); Jx_mass = reshape(Jx_mass,[total_time_step,1]);
    %Jy_mass = J_mass(:,2,:); Jy_mass = reshape(Jy_mass,[total_time_step,1]);
    %Jz_mass = J_mass(:,3,:); Jz_mass = reshape(Jz_mass,[total_time_step,1]);
    
    %Jx_all = Jx_all(1:total_time_step+1);
    %Jy_all = Jy_all(1:total_time_step+1);
    %Jz_all = Jz_all(1:total_time_step+1); 
    %Jcx_all= Jcx_all(1:total_time_step+1); 
    %Jcy_all= Jcy_all(1:total_time_step+1); 
    %Jcz_all= Jcz_all(1:total_time_step+1);
    Ns = Time_all(2) - Time_all(1); % sample interval in lammps

    %Jx_all = Jpx_all + Jcx_all; Jy_all = Jpy_all + Jcy_all; Jz_all = Jpz_all + Jcz_all;
    Jpx_all = Jx_all - Jcx_all; Jpy_all = Jy_all - Jcy_all; Jpz_all = Jz_all - Jcz_all;
    %Jz_all = Jz_all - mean(Jz_all);
    %Jpz_all = Jpz_all - mean(Jpz_all);
    %Jcz_all = Jcz_all - mean(Jcz_all);
    %Jz_mass = Jz_mass - mean(Jz_mass);

    
    finish = 'Begin to calculate hacf';
    % calculate heat current correlation function
    HAC_all = correc(Nc_all,Jx_all,Jy_all,Jz_all,Jx_all,Jy_all,Jz_all);
    HAC_ct = correc(Nc_all,Jcx_all,Jcy_all,Jcz_all,Jcx_all,Jcy_all,Jcz_all);
    HAC_pt = correc(Nc_all,Jpx_all,Jpy_all,Jpz_all,Jpx_all,Jpy_all,Jpz_all);
    HAC_cpt = correc(Nc_all,Jcx_all,Jcy_all,Jcz_all,Jpx_all,Jpy_all,Jpz_all);
    HAC_pct = correc(Nc_all,Jpx_all,Jpy_all,Jpz_all,Jcx_all,Jcy_all,Jcz_all);

    hac_all = [hac_all,HAC_all];hac_ct = [hac_ct,HAC_ct];hac_pt = [hac_pt,HAC_pt];
    hac_cpt = [hac_cpt,HAC_cpt];hac_pct = [hac_pct,HAC_pct];

    % calculate mass current correlation function
    m_Ag = 107;
    m_Te = 127;
    %number_Ag = 1728;
    h_Ag_mean = mean(h_Ag);
    h_Te_mean = mean(h_Te);
    h_Ag_Te = h_Ag_mean/m_Ag-h_Te_mean/m_Te;
    %h_Ag_Te = h_Ag./m_Ag-h_Te./m_Te;
    %vx = Jx_mass/m_Ag;
    %vy = Jy_mass/m_Ag;
    %vz = Jz_mass/m_Ag;

    Jx_mass_h = Jx_mass.*h_Ag_Te;
    Jy_mass_h = Jy_mass.*h_Ag_Te;
    Jz_mass_h = Jz_mass.*h_Ag_Te;
    Jx_cond = Jx_all - Jx_mass_h;
    Jy_cond = Jy_all - Jy_mass_h;
    Jz_cond = Jz_all - Jz_mass_h;

    HAC_cond = correc(Nc_all,Jx_cond,Jx_cond,Jx_cond,Jx_cond,Jx_cond,Jx_cond);
    hac_cond = [hac_cond,HAC_cond];

    MACE_all = correc(Nm_all,Jx_mass_h,Jy_mass_h,Jz_mass_h,Jx_mass_h,Jy_mass_h,Jz_mass_h);
    MAC_all = correc(Nm_all,Jx_mass,Jy_mass,Jz_mass,Jx_mass,Jy_mass,Jz_mass);
    mac_all = [mac_all,MAC_all];
    mace_all = [mace_all,MACE_all];

    % calculate mass/heat current correlation function
    HMACE_cond = correc(Ncm_all,Jx_mass_h,Jy_mass_h,Jz_mass_h,Jx_cond,Jx_cond,Jx_cond);
    HMACE_cond_2 = correc(Ncm_all,Jx_cond,Jx_cond,Jx_cond,Jx_mass_h,Jy_mass_h,Jz_mass_h);
    HMAC_all = correc(Ncm_all,Jx_mass,Jy_mass,Jz_mass,Jx_all,Jy_all,Jz_all);
    HMAC_cond = correc(Ncm_all,Jx_mass,Jy_mass,Jz_mass,Jx_cond,Jy_cond,Jz_cond);
    hmac_all = [hmac_all,HMAC_all];
    hmac_cond = [hmac_cond,HMAC_cond];
    hmace_cond = [hmace_cond,HMACE_cond];
    hmace_cond_2 = [hmace_cond_2,HMACE_cond_2];

    %Diffusion_coe = correc(Nm_all,vx,vy,vz,vx,vy,vz);
    %diffusion_coe = [diffusion_coe,Diffusion_coe];

    %HAC_fft = fft(HAC_all);
    %hac_fft = [hac_fft,HAC_fft];
    %MAC_fft  = fft(MAC_all);
    %mac_fft = [mac_fft,MAC_fft];
    %HMAC_fft = fft(HMAC_all);
    %hmac_fft = [hmac_fft,HMAC_fft];
end
    
finish = 'Begin to compute RTC';
% calculate thermal conductivity
dt = time_step*Ns;                                            % time for intergratioin
T = mean(Temperature(1:end-1, :));                            % temperature in every simulation
scale_c = 1./(kB .* T.* T.* V)* dt*ev2J/ps2s/a2m;                             % scale for unit change
scale_m = 1./(kB .* T.* T.* V)* dt/ev2J/ps2s/a2m*g2kg*g2kg/NA/NA;
scale_cm = 1./(kB .* T.* T.* V)* dt/ps2s/a2m*g2kg/NA;
%% correlation time
tc = (0 : Nc_all-1) * dt; % correlation time
tm = (0 : Nm_all-1) * dt;
tcm =(0 : Ncm_all-1) * dt;


%calculate RTC from hac

%%% total 
rtc_all =  scale_c .* cumtrapz(hac_all);
rtc_ct =  scale_c .* cumtrapz(hac_ct);
rtc_pt =  scale_c .* cumtrapz(hac_pt);
rtc_cpt =  scale_c .* cumtrapz(hac_cpt);
rtc_pct =  scale_c .* cumtrapz(hac_pct);
rtc_cond = scale_c .* cumtrapz(hac_cond);
rtme = scale_c .* cumtrapz(mace_all);
rtcme = scale_c .* cumtrapz(hmace_cond);
rtcme_2 = scale_c .* cumtrapz(hmace_cond_2);

rtm = scale_m .* cumtrapz(mac_all);
rtcm = scale_cm .* cumtrapz(hmac_all);
rtcm_cond = scale_cm .* cumtrapz(hmac_cond);

%hac_fft = scale_c .*mean(hac_fft,2);
%mac_fft = scale_m .*mean(mac_fft,2);
%hmac_fft = scale_cm .*mean(hmac_fft,2);

%%% average and acquire the final result
rtc_all_ave = mean(rtc_all, 2); 
rtc_ct_ave = mean(rtc_ct, 2);
rtc_pt_ave = mean(rtc_pt, 2);
rtc_cpt_ave = mean(rtc_cpt, 2);
rtc_pct_ave = mean(rtc_pct, 2);

rtc_cond_ave = mean(rtc_cond, 2);
rtm_ave = mean(rtm, 2);
rtcm_ave = mean(rtcm, 2);
rtcm_cond_ave = mean(rtcm_cond, 2);

rtme_ave = mean(rtme,2);
rtcme_ave = mean(rtcme,2);
rtcme_ave_2 = mean(rtcme_2,2);
%diffusion_coe_t = mean(cumtrapz(diffusion_coe),2);

LEE = mean(rtc_all_ave(Nc_all/10*6+1:Nc_all));
LE1 = mean(rtcm_ave(Nc_all/10*6+1:Nc_all));
L11 = mean(rtm_ave(Nc_all/10*6+1:Nc_all));
%result_E = [LEE LE1*LE1/L11 LEE-LE1*LE1/L11];

k_c = mean(rtc_ct_ave(Nc_all/10*6+1:Nc_all));
k_p = mean(rtc_pt_ave(Nc_all/10*6+1:Nc_all));
k_cp = mean(rtc_cpt_ave(Nc_all/10*6+1:Nc_all));
k_pc = mean(rtc_pct_ave(Nc_all/10*6+1:Nc_all));

k_c_error = std(mean(rtc_ct(Nc_all/10*6+1:Nc_all,:),1)) / sqrt(NN);
k_p_error = std(mean(rtc_pt(Nc_all/10*6+1:Nc_all,:),1)) / sqrt(NN);
k_cp_error = std(mean(rtc_cpt(Nc_all/10*6+1:Nc_all,:),1)) / sqrt(NN);
k_pc_error = std(mean(rtc_pct(Nc_all/10*6+1:Nc_all,:),1)) / sqrt(NN);

Lqq = mean(rtc_cond_ave(Nc_all/10*6+1:Nc_all));
Lq1 = mean(rtcm_cond_ave(Ncm_all/10*6+1:Ncm_all));
result_q = [Lqq Lq1*Lq1/L11 Lqq-Lq1*Lq1/L11];

L11h = mean(rtme_ave(Nc_all/10*6+1:Nc_all));
Lq1h = mean(rtcme_ave(Nc_all/10*6+1:Nc_all));
L1qh = mean(rtcme_ave_2(Nc_all/10*6+1:Nc_all));

LEE_error = std(mean(rtc_all(Nc_all/10*6+1:Nc_all,:),1)) / sqrt(NN);
Lqq_error = std(mean(rtc_cond(Nc_all/10*6+1:Nc_all,:),1)) / sqrt(NN);
Lq1h_error = std(mean(rtcme(Nc_all/10*6+1:Nc_all,:),1)) / sqrt(NN);
L1qh_error = std(mean(rtcme_2(Nc_all/10*6+1:Nc_all,:),1)) / sqrt(NN);
L11h_error = std(mean(rtme(Nc_all/10*6+1:Nc_all,:),1)) / sqrt(NN);

%result_E = [LEE LEE_error LE1*1e6 L11*1e12 LEE-LE1*LE1/L11]
%result_qh = [Lqq Lqq_error L11h L11h_error Lq1h Lq1h_error Lqq-Lq1*Lq1/L11]
reduced_unit = 0.1250;
result_qh = [LEE LEE_error Lqq Lqq_error L11h L11h_error Lq1h+L1qh Lq1h_error+L1qh_error]./reduced_unit
result_k = [LEE LEE_error k_p k_p_error k_c k_c_error k_cp+k_pc k_pc_error+k_cp_error]./reduced_unit
%diffusion_coe_mean = mean(diffusion_coe_t(Nc_all/10*6+1:Nc_all))

%% unit transfer


%kappa_all_ke_error = std(mean(rtc_all_ke(Nc_all/10*6+1:Nc_all,:),1)) / sqrt(NN);
%kappa_all_pe_error = std(mean(rtc_all_pe(Nc_all/10*6+1:Nc_all,:),1)) / sqrt(NN);

%result_Lqq = [Lqq_all Lqq_c Lqq_p Lqq_cp Lqq_pc (Lqq_cp+Lqq_pc)] 
%result_L11 = L11
%result_L11_H = L11*H_Ag_Te*H_Ag_Te
%result_L1q = L1q
%result_L1q_H = L1q*H_Ag_Te*2
%result_L11_L1q = L1q*L1q/L11

tc = tc';
tm = tm';
tcm = tcm';

HAC_all = HAC_all/HAC_all(1);
HAC_cond = HAC_cond/HAC_cond(1);
MACE_all = MACE_all/MACE_all(1);

subplot(2,1,1)
plot(tc,HAC_all,'black--')
hold on
plot(tc,HAC_cond,'red--')
hold on
plot(tm,MACE_all,'blue--')
%hold on
%plot(tcm,HMACE_cond,'black--')

%subplot(2,1,2)
%plot(tc,diffusion_coe_t,'red--')
subplot(2,1,2)
plot(tc,rtc_cond,'black--')
hold on
plot(tc,rtme,'red--')
hold on
plot(tm,rtcme,'blue--')

function HAC = correc(Nc,J1x,J1y,J1z,J2x,J2y,J2z)
Nd = size(J1x, 1);
m = Nd - Nc;
HACx = [];HACy = [];HACz = [];
for i = 0:Nc-1
   ax = 0; ay = 0; az = 0;   % <Jx*Jx> ... the first term is Jx(t')
   for j = 1:m
       ax = ax + J1x(j)*J2x(i+j);   ay = ay + J1y(j)*J2y(i+j);   az = az + J1z(j)*J2z(i+j);
   end
   HACx = [HACx; ax/m];     HACy = [HACy; ay/m];     HACz = [HACz; az/m];
end
HAC = (HACx + HACy + HACz) / 3;
%HAC = HACz;
end
