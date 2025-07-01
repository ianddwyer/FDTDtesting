
%holdOutput is for each of the ports which each have
%sub index of the freq and time for all possible probes
%use repetition to pull all values and use line test as ref values

obj_180hybrid = struct([]);
field_data = struct([]);
netparams = struct([]);
linespecs = struct([]);
mesh(1).dt = dt;
nopFreq = floor(opFreq/df);
try
    close all;
    load("line_wkspc.mat") %%%MUST HAVE 180HYBRID SAVED IN POSTPROCESS TO PREVENT LOSING DATA BEFORE LOAD!!!
    config_line; %regenerate PCB as figure 2 (hardcoded in hybrid)
    obj_180hybrid(1).object = "line";
    obj_180hybrid(1).outputs = holdOutput_line; %load line object holds after line runs in try...
    obj_180hybrid(1).workspace = "line_wkspc.mat"; %output data from time and freq domain
    obj_180hybrid(1).PEC = PEC; %the trace's coordinates from orientation object definitions
    obj_180hybrid(1).trace = trace; %orientation of trace
    obj_180hybrid(1).linespecs = linespecs; %lengths of trace relations
    obj_180hybrid(1).mesh = mesh; %environment, physical and time domain setup
    obj_180hybrid(1).source = source; %source data
    obj_180hybrid(1).port = port; %port orientation data
    obj_180hybrid(1).probe = probe; %probe reference information
    obj_180hybrid(1).substrate = substrate; %the dielectric used for testing
    obj_180hybrid(1).opFreq = opFreq; %operational frerquency of test
    obj_180hybrid(1).quarterwav_mil = qwav; %quarter wavelength of circuit in mil

    load 180hybrid_wkspc.mat
    config_180hybrid; %regenerate PCB as figure 1 (hardcoded in line) 
catch
    %
end
%load ALL of to ensure correct one is loaded first due to same var names
Vhybrid = holdOutput(1).time(1).data;%from port 1 of hybrid to compare to port 1 at source on 1, hold before load
%only load certain values from line for reference
load line_wkspc.mat Z holdOutput_line OutputValue; 
Zo = Z;
%known forward wave from line test

figure
set(gcf,"name",'NETWORK PARAMETERS')
set(gcf,'Position',[100+floor(900*3/4) 100 floor(900*3/4) floor(720*3/4)])

opFreq=3e9;

for ss=sourceport %sourceport is array and sourcePort is scalar
    
    Vline = holdOutput_line(1).time(1).data; %line probe order is V I I
    %subtract forward wave from total wave to obtain reverse wave
    Vrev = Vhybrid;
    Vrev(98:380) = (Vhybrid(98:380)-Vline(98:380));%accelerated code requires scaling
    holdOutput(ss).freq(1).data = 20.*log10(fft(holdOutput(ss).time(1).data,NFFT)./scalingV_coeff);
    holdOutput(ss).freq(2).data = 20.*log10(fft(holdOutput(ss).time(2).data,NFFT)./scalingV_coeff);
    holdOutput(ss).freq(3).data = 20.*log10(fft(holdOutput(ss).time(3).data,NFFT)./scalingV_coeff);
    holdOutput(ss).freq(4).data = 20.*log10(fft(holdOutput(ss).time(4).data,NFFT)./scalingV_coeff);
    holdOutput(ss).freq(ss).data = 20.*log10(fft(Vrev,NFFT)./((scalingV_coeff)));

    subplot(220+ss)

    plot(holdOutput(ss).f,real(holdOutput(ss).freq(1).data)); 
    hold on
    plot(holdOutput(ss).f,real(holdOutput(ss).freq(2).data))
    plot(holdOutput(ss).f,real(holdOutput(ss).freq(3).data))
    plot(holdOutput(ss).f,real(holdOutput(ss).freq(4).data))
    title("S Network Params w/ Gaus Source at Port "+string(ss))
    ylabel("(V/V)^2 [dB]")
    xlabel("frequency[Hz]")
    Straces = "S"+["1","2","3","4"]+string(ss);
    legend(Straces, Location="southeast")
    xlim([0 6]*1e9)
    ylim([-50 3])
    grid on

end





SFDTD = [...
    holdOutput(1).freq(1).data(nopFreq),holdOutput(1).freq(2).data(nopFreq),holdOutput(1).freq(3).data(nopFreq),holdOutput(1).freq(4).data(nopFreq);...
    holdOutput(2).freq(1).data(nopFreq),holdOutput(2).freq(2).data(nopFreq),holdOutput(2).freq(3).data(nopFreq),holdOutput(2).freq(4).data(nopFreq);...
    holdOutput(3).freq(1).data(nopFreq),holdOutput(3).freq(2).data(nopFreq),holdOutput(3).freq(3).data(nopFreq),holdOutput(3).freq(4).data(nopFreq);...
    holdOutput(4).freq(1).data(nopFreq),holdOutput(4).freq(2).data(nopFreq),holdOutput(4).freq(3).data(nopFreq),holdOutput(4).freq(4).data(nopFreq)...
    ];
S1 = [...
    holdOutput(1).freq(1).data(:)';...
    holdOutput(2).freq(1).data(:)';...
    holdOutput(3).freq(1).data(:)';...
    holdOutput(4).freq(1).data(:)';...
    ];
S2 = [...
    holdOutput(1).freq(2).data(:)';...
    holdOutput(2).freq(2).data(:)';...
    holdOutput(3).freq(2).data(:)';...
    holdOutput(4).freq(2).data(:)';...
    ];
S3 = [...
    holdOutput(1).freq(3).data(:)';...
    holdOutput(2).freq(3).data(:)';...
    holdOutput(3).freq(3).data(:)';...
    holdOutput(4).freq(3).data(:)';...
    ];
S4 = [...
    holdOutput(1).freq(4).data(:)';...
    holdOutput(2).freq(4).data(:)';...
    holdOutput(3).freq(4).data(:)';...
    holdOutput(4).freq(4).data(:)';...
    ];

Sij = zeros([4,4,length(holdOutput(1).freq(2).data(:))]);
Sij(:,1,:)=S1;
Sij(:,2,:)=S2;
Sij(:,3,:)=S3;
Sij(:,4,:)=S4;
% holdSigns = sign(imag(SFDTD)).*sign(real(SFDTD)); %due to finding Sii imag will always be 0, add identity matrix to prevent zeroing when reapplying the signs
% holdSigns = -holdSigns; %signs flip and real imaginary trade places due to being multiplied by j

SFDTD = abs( 10.^(SFDTD/20)*sqrt(2)*1j );%.*holdSigns; %normalized S-parameters
%     error = ([ 0 -1 -1 0; -1 0 0 1; -1 0 0 -1; 0 1 -1 0]-S)./[ 0 -1 -1 0; -1 0 0 1; -1 0 0 -1; 0 1 -1 0]


%% Other Params

%% effective dielec from impedence to optimize PML

try
    close([10,11])
catch
    %
end

% if config=="line"
    beta_eff = -angle( ( fft(OutputValue(:,1),NFFT)./fft(OutputValue(:,2),NFFT))  )./((shift_probe_from_source)*2*dy);
    omega = 2*pi.*holdOutput_line(1).f; 
    eps_eff = ((beta_eff).^2./(omega*c0^-1)'.^2);%./(mesh(1).dx/2); %this scaling is wrong... completely but value seems right....
    % disp(eps_eff(round(opFreq/df)))

    figure(10); 
    plot(holdOutput_line(1).f,eps_eff)
    xlim([0 6e9])
    ylim([0 12])
    grid on
    title("Effective Substrate Dielectric Permitivity")
    xlabel("frequency [Hz]")
    ylabel("\epsilon_r [F/m]")


%     figure(11)
% %     plot(OutputValue(92:385,1))
%     hold on
% %     plot(OutputValue(92:385,3))
%     hold on
%     Jy = OutputValue(92:385,4)./(dy/2)./(Zo(round(opFreq./df)));
%     Vz = Jy*dt/(epsrB*eps0)/(dz*subheight)./4;
%     plot(holdOutput_line(1).t,Vz)
%     legend("Vz")
%     grid on

% end
%%

Z0 = holdOutput_line(1).Z(nopFreq)
ABCD = s2abcd(SFDTD,Z0);
Zmat = abs(s2z(-abs(SFDTD*(1/sqrt(2))),Z0))
Ymat = pageinv(Zmat);


netparams(1).eps_eff = eps_eff(round(opFreq./df));
netparams(1).Smatrix = SFDTD;
netparams(1).Ymatrix = Ymat;
netparams(1).Zmatrix = Zmat;
netparams(1).ABCD = ABCD;


field_data(1).x = ex;
field_data(1).y = ey;
field_data(1).z = ez;
field_data(2).x = hx;
field_data(2).y = ey;
field_data(2).z = ez;

%object 1 is line data, object 2 is 180hybrid data
obj_180hybrid(2).object = config; %object type for reference
obj_180hybrid(2).workspace = "180hybrid_wkspc.mat"; %output data from time and freq domain
obj_180hybrid(2).outputs = holdOutput; %output data from time and freq domain
obj_180hybrid(2).netparams = netparams; %scattering paramters
obj_180hybrid(2).fields = field_data; %output of final time field data
obj_180hybrid(2).PEC = PEC; %the trace's coordinates from orientation object definitions
obj_180hybrid(2).trace = trace; %orientation of trace
obj_180hybrid(2).linespecs = linespecs; %lengths of trace relations
obj_180hybrid(2).mesh = mesh; %environment, physical and time domain setup
obj_180hybrid(2).source = source; %source data
obj_180hybrid(2).port = port; %port orientation data
obj_180hybrid(2).probe = probe; %probe reference information
obj_180hybrid(2).substrate = substrate; %the dielectric used for testing
obj_180hybrid(2).opFreq = opFreq; %operational frerquency of test
obj_180hybrid(2).quarterwav_mil = qwav; %quarter wavelength of circuit in mil

%% ADS data
% S_ADS = sparameters('hybrid_data_ADS_MoM.s4p');
S_ADS = sparameters('hybrid_data_ADS_spice.s4p');
% S_ADS = sparameters('MoMResults.s4p');
NfMax = floor((6e9)/df);
fNew = (0:NfMax-1)*df;
fADS = S_ADS(1).Frequencies;
datADS = squeeze(S_ADS(1).Parameters(1,1,:));
S11_spice = 20*log10(abs(interp1(fADS,datADS,fNew)));
datADS = squeeze(S_ADS(1).Parameters(1,2,:));
S12_spice = 20*log10(abs(interp1(fADS,datADS,fNew)));
datADS = squeeze(S_ADS(1).Parameters(1,3,:));
S13_spice = 20*log10(abs(interp1(fADS,datADS,fNew)));
datADS = squeeze(S_ADS(1).Parameters(1,4,:));
S14_spice = 20*log10(abs(interp1(fADS,datADS,fNew)));
S11_FDTD=real(holdOutput(1).freq(1).data(1:NfMax));
S12_FDTD=real(holdOutput(2).freq(1).data(1:NfMax));
S13_FDTD=real(holdOutput(3).freq(1).data(1:NfMax));
S14_FDTD=real(holdOutput(4).freq(1).data(1:NfMax));
try
    close(7)
catch
    %
end
figure(7); 
subplot(221)
plot(fNew,S11_spice)
hold on
plot(fNew,S12_spice)
plot(fNew,S13_spice)
plot(fNew,S14_spice)
plot(fNew,S11_FDTD,':',LineWidth=2); 
plot(fNew,S12_FDTD,':',LineWidth=2);
plot(fNew,S13_FDTD,':',LineWidth=2);
plot(fNew,S14_FDTD,':',LineWidth=2);
title("Spice (ADS) vs FDTD (Custom) for 180-hybrid Coupler")
xlim([0 6]*1e9)
ylim([-70 5])
legend("spice S11","spice S12","spice S13","spice S14", "FDTD S11","FDTD S12","FDTD S13","FDTD S14", Location="southeast")
yline(0)
grid on
ylabel("Port Gain [dB]")
xlabel("Frequency [Hz]")

subplot(223)
diffS11 = abs((  S11_FDTD - S11_spice'  ));
diffS12 = abs((  S12_FDTD - S12_spice'  ));
diffS13 = abs((  S13_FDTD - S13_spice'  ));
diffS14 = abs((  S14_FDTD - S14_spice'  ));
plot(fNew,diffS11)
hold on
plot(fNew,diffS12)
plot(fNew,diffS13)
plot(fNew,diffS14)
% scatter(S11_FDTD',S11_spice)
ylim([0 70])
title("Absolute Difference: Spice vs FDTD")
legend("S11","S12","S13","S14")
ylabel("Diff [dB]")
xlabel("frequency [Hz]")
grid on


%S_ADS = sparameters('hybrid_data_ADS_MoM.s4p');
%S_ADS = sparameters('hybrid_data_ADS_spice.s4p');
S_ADS = sparameters('MoMResults.s4p');
NfMax = floor((6e9)/df);
fNew = (0:NfMax-1)*df;
fADS = S_ADS(1).Frequencies;
datADS = squeeze(S_ADS(1).Parameters(1,1,:));
S11_MoM = 20*log10(abs(interp1(fADS,datADS,fNew)));
datADS = squeeze(S_ADS(1).Parameters(1,2,:));
S12_MoM = 20*log10(abs(interp1(fADS,datADS,fNew)));
datADS = squeeze(S_ADS(1).Parameters(1,3,:));
S13_MoM = 20*log10(abs(interp1(fADS,datADS,fNew)));
datADS = squeeze(S_ADS(1).Parameters(1,4,:));
S14_MoM = 20*log10(abs(interp1(fADS,datADS,fNew)));
S11_FDTD=real(holdOutput(1).freq(1).data(1:NfMax));
S12_FDTD=real(holdOutput(1).freq(2).data(1:NfMax));
S13_FDTD=real(holdOutput(1).freq(3).data(1:NfMax));
S14_FDTD=real(holdOutput(1).freq(4).data(1:NfMax));

subplot(222)
plot(fNew,S11_MoM)
hold on
plot(fNew,S12_MoM)
plot(fNew,S13_MoM)
plot(fNew,S14_MoM)
plot(fNew,S11_FDTD,':',LineWidth=2); 
plot(fNew,S12_FDTD,':',LineWidth=2);
plot(fNew,S13_FDTD,':',LineWidth=2);
plot(fNew,S14_FDTD,':',LineWidth=2);
title("MoM (ADS) vs FDTD (Custom) for 180-hybrid Coupler")
xlim([0 6]*1e9)
ylim([-70 5])
legend("MoM S11","MoM S12","MoM S13","MoM S14", "FDTD S11","FDTD S12","FDTD S13","FDTD S14", Location="southeast")
yline(0)
grid on
ylabel("Port Gain [dB]")
xlabel("Frequency [Hz]")

subplot(224)
diffS11 = abs((  S11_FDTD - S11_MoM'  ));
diffS12 = abs((  S12_FDTD - S12_MoM'  ));
diffS13 = abs((  S13_FDTD - S13_MoM'  ));
diffS14 = abs((  S14_FDTD - S14_MoM'  ));
plot(fNew,diffS11)
hold on
plot(fNew,diffS12)
plot(fNew,diffS13)
plot(fNew,diffS14)
% scatter(S11_FDTD',S11_MoM)
ylim([0 70])
title("Absolute Difference: MoM vs FDTD")
legend("S11","S12","S13","S14")
ylabel("Diff [dB]")
xlabel("frequency [Hz]")
grid on








% subplot(313)
% diffS11 = abs((  S11_FDTD - S11_MoM'  )./S11_FDTD)*100;
% diffS12 = abs((  S12_FDTD - S12_MoM'  )./S11_FDTD)*100;
% diffS13 = abs((  S13_FDTD - S13_MoM'  )./S11_FDTD)*100;
% diffS14 = abs((  S14_FDTD - S14_MoM'  )./S11_FDTD)*100;
% plot(fNew,diffS11)
% hold on
% plot(fNew,diffS12)
% plot(fNew,diffS13)
% plot(fNew,diffS14)
% % scatter(S11_FDTD',S11_MoM)
% title("Relative Difference: MoM vs FDTD")
% legend("S11","S12","S13","S14")
% ylabel("Diff [%]")
% xlabel("frequency [Hz]")
% ylim([0 100])
% grid on

% 
% clr_wkspc;
% 
% subplot(312)
% plot(holdOutput(1).f,holdOutput_line(1).Z)
% xlim([0 6]*1e9)
% grid on
% legend("Z")
% subplot(313)
% plot(holdOutput(1).f,gamma)
% legend("gamma")
% yline(0)
% grid on
% xlim([0 6]*1e9)

%% EDGE LENGTH AND H_FIELD SCALING TO FIX THESE ERRORS!!!!
%% a optimal calculations in notes to optimize PML

%% phase
figure
set(gcf,"name",'NETWORK PARAMETERS')
set(gcf,'Position',[100+floor(900*3/4) 100 floor(900*3/4) floor(720*3/4)])

opFreq=3e9;

for ss=sourceport %sourceport is array and sourcePort is scalar
    
    Vline = holdOutput_line(1).time(1).data; %line probe order is V I I
    %subtract forward wave from total wave to obtain reverse wave
    Vrev = (Vhybrid-Vline);%accelerated code requires scaling
    holdOutput(ss).freq(1).data = abs(fft(holdOutput(ss).time(1).data,NFFT)./scalingV_coeff);
    holdOutput(ss).freq(2).data = abs(fft(holdOutput(ss).time(2).data,NFFT)./scalingV_coeff);
    holdOutput(ss).freq(3).data = abs(fft(holdOutput(ss).time(3).data,NFFT)./scalingV_coeff);
    holdOutput(ss).freq(4).data = abs(fft(holdOutput(ss).time(4).data,NFFT)./scalingV_coeff);
    holdOutput(ss).freq(ss).data = abs(fft(Vrev,NFFT)./((scalingV_coeff)));

    subplot(220+ss)

    plot(holdOutput(ss).f,abs(holdOutput(ss).freq(1).data)); 
    hold on
    plot(holdOutput(ss).f,abs(holdOutput(ss).freq(2).data))
    plot(holdOutput(ss).f,abs(holdOutput(ss).freq(3).data))
    plot(holdOutput(ss).f,abs(holdOutput(ss).freq(4).data))
    title("S Network Params w/ Gaus Source at Port "+string(ss))
    ylabel("(V/V)^2 [rad]")
    xlabel("frequency[Hz]")
    Straces = "S"+["1","2","3","4"]+string(ss);
    legend(Straces, Location="southeast")
    xlim([0 6]*1e9)
    ylim([0 1])
    grid on

end

%% ADS data
% S_ADS = sparameters('hybrid_data_ADS_MoM.s4p');
S_ADS = sparameters('hybrid_data_ADS_spice.s4p');
% S_ADS = sparameters('MoMResults.s4p');
NfMax = floor((6e9)/df);
fNew = (0:NfMax-1)*df;
fADS = S_ADS(1).Frequencies;
datADS = squeeze(S_ADS(1).Parameters(1,1,:));
S11_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(1,2,:));
S12_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(1,3,:));
S13_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(1,4,:));
S14_spice = abs(interp1(fADS,datADS,fNew));
S11_FDTD=real(holdOutput(1).freq(1).data(1:NfMax));
S12_FDTD=real(holdOutput(2).freq(1).data(1:NfMax));
S13_FDTD=real(holdOutput(3).freq(1).data(1:NfMax));
S14_FDTD=real(holdOutput(4).freq(1).data(1:NfMax));
try
    close(7)
catch
    %
end
figure(7); 
subplot(221)
plot(fNew,20.*log10(S11_spice))
hold on
plot(fNew,20.*log10(S12_spice))
plot(fNew,20.*log10(S13_spice))
plot(fNew,20.*log10(S14_spice))
plot(fNew,20.*log10(S11_FDTD),':',LineWidth=2); 
plot(fNew,20.*log10(S12_FDTD),':',LineWidth=2);
plot(fNew,20.*log10(S13_FDTD),':',LineWidth=2);
plot(fNew,20.*log10(S14_FDTD),':',LineWidth=2);
title("Spice (ADS) vs FDTD (Custom) for 180-hybrid Coupler")
xlim([0 6]*1e9)
% ylim([-70 5])
legend("spice S11","spice S12","spice S13","spice S14", "FDTD S11","FDTD S12","FDTD S13","FDTD S14", Location="northwest")
yline(0)
grid on
ylabel("Port Gain [dB]")
xlabel("Frequency [Hz]")

subplot(223)
diffS11 = abs((  S11_FDTD - S11_spice'  )).*100;
diffS12 = abs((  S12_FDTD - S12_spice'  )).*100;
diffS13 = abs((  S13_FDTD - S13_spice'  )).*100;
diffS14 = abs((  S14_FDTD - S14_spice'  )).*100;
plot(fNew,diffS11)
hold on
plot(fNew,diffS12)
plot(fNew,diffS13)
plot(fNew,diffS14)
% scatter(S11_FDTD',S11_spice)
ylim([0 100])
title("Absolute Difference: Spice vs FDTD")
legend("S11","S12","S13","S14",Location="northwest")
ylabel("Diff [%]")
xlabel("frequency [Hz]")
grid on


%S_ADS = sparameters('hybrid_data_ADS_MoM.s4p');
%S_ADS = sparameters('hybrid_data_ADS_spice.s4p');
S_ADS = sparameters('MoMResults.s4p');

NfMax = floor((6e9)/df);
fNew = (0:NfMax-1)*df;
fADS = S_ADS(1).Frequencies;
datADS = squeeze(S_ADS(1).Parameters(1,1,:));
S11_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(1,2,:));
S12_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(1,3,:));
S13_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(1,4,:));
S14_MoM = abs(interp1(fADS,datADS,fNew));
S11_FDTD=real(holdOutput(1).freq(1).data(1:NfMax));
S12_FDTD=real(holdOutput(1).freq(2).data(1:NfMax));
S13_FDTD=real(holdOutput(1).freq(3).data(1:NfMax));
S14_FDTD=real(holdOutput(1).freq(4).data(1:NfMax));

subplot(222)
plot(fNew,20.*log10(S11_MoM))
hold on
plot(fNew,20.*log10(S12_MoM))
plot(fNew,20.*log10(S13_MoM))
plot(fNew,20.*log10(S14_MoM))
plot(fNew,20.*log10(S11_FDTD),':',LineWidth=2); 
plot(fNew,20.*log10(S12_FDTD),':',LineWidth=2);
plot(fNew,20.*log10(S13_FDTD),':',LineWidth=2);
plot(fNew,20.*log10(S14_FDTD),':',LineWidth=2);
title("MoM (ADS) vs FDTD (Custom) for 180-hybrid Coupler")
xlim([0 6]*1e9)
% ylim([-70 5])
legend("MoM S11","MoM S12","MoM S13","MoM S14", "FDTD S11","FDTD S12","FDTD S13","FDTD S14", Location="northwest")
yline(0)
grid on
ylabel("Port Gain [dB]")
xlabel("Frequency [Hz]")

subplot(224)
diffS11 = abs((  S11_FDTD - S11_MoM'  )).*100;
diffS12 = abs((  S12_FDTD - S12_MoM'  )).*100;
diffS13 = abs((  S13_FDTD - S13_MoM'  )).*100;
diffS14 = abs((  S14_FDTD - S14_MoM'  )).*100;
plot(fNew,diffS11)
hold on
plot(fNew,diffS12)
plot(fNew,diffS13)
plot(fNew,diffS14)
% scatter(S11_FDTD',S11_MoM)
ylim([0 100])
title("Absolute Difference: MoM vs FDTD")
legend("S11","S12","S13","S14", Location="northwest")
ylabel("Diff [%]")
xlabel("frequency [Hz]")
grid on

epsEff = [eps_eff(round(opFreq./df))]

datADS = squeeze(S_ADS(1).Parameters(1,1,:));
S11_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(1,2,:));
S12_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(1,3,:));
S13_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(1,4,:));
S14_MoM = abs(interp1(fADS,datADS,fNew));

datADS = squeeze(S_ADS(1).Parameters(2,1,:));
S21_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(2,2,:));
S22_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(2,3,:));
S23_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(2,4,:));
S24_MoM = abs(interp1(fADS,datADS,fNew));

datADS = squeeze(S_ADS(1).Parameters(3,1,:));
S31_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(3,2,:));
S32_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(3,3,:));
S33_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(3,4,:));
S34_MoM = abs(interp1(fADS,datADS,fNew));

datADS = squeeze(S_ADS(1).Parameters(4,1,:));
S41_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(4,2,:));
S42_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(4,3,:));
S43_MoM = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(4,4,:));
S44_MoM = abs(interp1(fADS,datADS,fNew));

% S_ADS = sparameters('hybrid_data_ADS_MoM.s4p');
S_ADS = sparameters('hybrid_data_ADS_spice.s4p');
% S_ADS = sparameters('MoMResults.s4p');

NfMax = floor((6e9)/df);
fNew = (0:NfMax-1)*df;
fADS = S_ADS(1).Frequencies;
datADS = squeeze(S_ADS(1).Parameters(1,1,:));
S11_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(1,2,:));
S12_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(1,3,:));
S13_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(1,4,:));
S14_spice = abs(interp1(fADS,datADS,fNew));

datADS = squeeze(S_ADS(1).Parameters(2,1,:));
S21_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(2,2,:));
S22_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(2,3,:));
S23_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(2,4,:));
S24_spice = abs(interp1(fADS,datADS,fNew));

datADS = squeeze(S_ADS(1).Parameters(3,1,:));
S31_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(3,2,:));
S32_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(3,3,:));
S33_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(3,4,:));
S34_spice = abs(interp1(fADS,datADS,fNew));

datADS = squeeze(S_ADS(1).Parameters(4,1,:));
S41_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(4,2,:));
S42_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(4,3,:));
S43_spice = abs(interp1(fADS,datADS,fNew));
datADS = squeeze(S_ADS(1).Parameters(4,4,:));
S44_spice = abs(interp1(fADS,datADS,fNew));

SFDTD = SFDTD
SMoM = ([...
    S11_MoM(nopFreq),S12_MoM(nopFreq),S13_MoM(nopFreq),S14_MoM(nopFreq);...
    S21_MoM(nopFreq),S22_MoM(nopFreq),S23_MoM(nopFreq),S24_MoM(nopFreq);...
    S31_MoM(nopFreq),S32_MoM(nopFreq),S33_MoM(nopFreq),S34_MoM(nopFreq);...
    S41_MoM(nopFreq),S42_MoM(nopFreq),S43_MoM(nopFreq),S44_MoM(nopFreq)...
    ].*sqrt(2))


Sspice = imag([...
    S11_spice(nopFreq),S12_spice(nopFreq),S13_spice(nopFreq),S14_spice(nopFreq);...
    S21_spice(nopFreq),S22_spice(nopFreq),S23_spice(nopFreq),S24_spice(nopFreq);...
    S31_spice(nopFreq),S32_spice(nopFreq),S33_spice(nopFreq),S34_spice(nopFreq);...
    S41_spice(nopFreq),S42_spice(nopFreq),S43_spice(nopFreq),S44_spice(nopFreq)...
    ].*sqrt(2).*1j)

%% BETTER METHOD TO RETRIEVE DATA BUT NOT ENOUGH TIME TO CLEAN UP PREV CODE
%PHASE DATA
try
close(12)
catch
end
S_ADS = sparameters('hybrid_data_ADS_MoM.s4p');
Sm = zeros([4 4 length(f)]);
clear temp
for row =1:4
    for col =1:4
        temp = squeeze(S_ADS(1).Parameters(row,col,:));
        Sm(row,col,:) = (interp1(fADS,temp,f));
    end
end
clear temp
S_ADS = sparameters('hybrid_data_ADS_MoM.s4p');
% S_ADS = sparameters('MoMResults.s4p');
Ss = zeros([4 4 length(f)]);
for row =1:4
    for col =1:4
        temp = squeeze(S_ADS(1).Parameters(row,col,:));
        Ss(row,col,:) = (interp1(fADS,temp,f));
    end
end

figure(12);
for jj = 1 %select source port
    for ii=1:4
        subplot(220+ii)
        if ii==jj 
            subvf = holdOutput_line(1).time(1).data; %remove forward signal
        else
            subvf = 0;
        end
        unwrapped_FDTD = 180/pi*unwrap(angle(fft(holdOutput(jj).time(ii).data-subvf,NFFT)./(NFFT*dx)));
        unwrapped_MoM = 180/pi*unwrap(angle(squeeze(Sm(jj,ii,:))));
        unwrapped_spice = 180/pi*unwrap(angle(squeeze(Ss(jj,ii,:))));
        phase_error_mean = sum(unwrapped_FDTD(1:nopFreq)-unwrapped_MoM(1:nopFreq))./length(1:nopFreq)
        phase_error_max = max(abs(unwrapped_FDTD(1:nopFreq)-unwrapped_MoM(1:nopFreq)))
        plot(holdOutput(1).f,unwrapped_FDTD,LineWidth=2);
        xlim([0 6e9]);
        % ylim([-360 0]);
        grid on; hold on; 
        
        mean_error = sum(unwrapped_FDTD(1:length(Ss(jj,ii,:))) - unwrapped_MoM)./(1:length(Ss(jj,ii,:)));
        plot(f,unwrapped_MoM,LineWidth=2)
        plot(f,unwrapped_spice,":",Color='k',LineWidth=1)
        legend("FDTD","MoM","Spice")
        xlabel("frequency [Hz]")
        ylabel('phase [deg]')%[\theta^\circ]')
        title("Phase Comparison: Source at Port "+string(jj)+" Load at Port "+string(ii))
    end
end