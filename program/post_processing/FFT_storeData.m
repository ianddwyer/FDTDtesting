%% fft for impedance and temporary output during design

NFFT = floor(2^(3*(log2(nMax)-9)+12));%calculates for cubic 
%df=pi/(2*dt*NFFT); %scales frequency for sample scaling in config
df=1/(dt*NFFT);
f = (0:NFFT-1).*df;
fftV = zeros([4 NFFT]);


% figure
Z=zeros([1 NFFT]);
Zall=zeros([length(sourceport) NFFT]);
Zo=zeros([1 NFFT]);
scalingV_coeff = NFFT*dx
scalingI_coeff = NFFT*dx
for probe_num = fftprobe %runs through requested port probes from app, suggested to use main for each port
    
    V1 = OutputValue(:,sourcePort);
    I1 = OutputValue(:,length(OutputValue(1,:))); %source curent probes are always last two probes!!
    I2 = OutputValue(:,length(OutputValue(1,:))-1); %source curent probes are always last two probes!!
    
    %gather source impedance                
    Z = abs(fft(V1./subheight,NFFT)./(sqrt(fft(I1,NFFT).*fft(I2,NFFT))));
    %gather impedance with respect to other ports
    count =1;
    if config~="line"  
        for ii = 1:2:length(sourceport)
            fV = fft(V1,NFFT)./subheight;
            fI1 = fft(OutputValue(:,ii+length(sourceport)),NFFT);
            fI2 = fft(OutputValue(:,ii+length(sourceport)+1),NFFT);
            Zall(count,:) = (  fV./sqrt(fI1.*fI2)  );
            count=count+1;
        end
    end

    if config~="line" %1 is stable...
        cd ..\..;
        load line_wkspc.mat Zo; 
        cd program\post_processing
        length(Zo)
    elseif config=="line"
        Zo = Z;
    end

    if sourcePort==probe(sourcePort).port || config=="line"
        parameter = "Port Impedance";
    else
        parameter = "Z_"+string(sourcePort)+"_"+string(probe(1).port);
    end

    sourceBlocks = subheight+1;%plus 1 due to unseen PEC ground from wall
    if probe(sourcePort).axis=="y"
        probeBlocks = (probe(sourcePort).j(2)-probe(1).j(1))+1;
    elseif probe(sourcePort).axis=="x"
        probeBlocks = (probe(sourcePort).i(2)-probe(1).i(1))+1;
    end
    
    fftV(probe_num,:)=fft(OutputValue(:,probe_num),NFFT)./scalingV_coeff;

    
%     plot(f,Z)
%     ylim([0 100])
%     xlim([0 10]*1e9)
%     xlabel('frequency [Hz]')
%     ylabel('Z(f)[\Omega]')
%     title("Frequency Magnitude Response for "+parameter)
%     grid on
    
    
%     plot(f,fftV(probe_num,:))
%     hold on
%     %plot(1/(pi*srct_w),fftV(floor(1/(pi*srct_w)/df)),'o') %corner of gaus
%     xlim([0 10].*1e9)
%     ylim([-80 80])
%     xlabel('frequency [Hz]')
%     ylabel('Power [dB]')
%     title("Frequency Magnitude Response with Port "+sourcePort+" Excited")
%     grid on
%     set(gcf,'units','normalized','OuterPosition',[0 0 1 1])
    
    fVecs(probe_num).data = fftV(probe_num,:);
    tVecs(probe_num).data = OutputValue(:,probe_num);
 
end

ports = zeros([1 nports+1]);
count = 1;
for ii = fftprobe
    ports(count) = probe(ii).port;
    count = count+1;
end

ports = string(ports);
ports = "port "+ports(1:nports);
ports(count) = "0dB";
yline(0)
legend(ports)


%store simulation data
if config ~= "line"
    holdOutput(sourcePort).time     = tVecs;
    holdOutput(sourcePort).freq     = fVecs;
    holdOutput(sourcePort).Z        = Z;
    holdOutput(sourcePort).Zports   = Zall;
    holdOutput(sourcePort).dt       = dt;
    holdOutput(sourcePort).Nt       = nMax;
    holdOutput(sourcePort).t        = dataOut(:,1);
    holdOutput(sourcePort).f        = f;
    holdOutput(sourcePort).Nf       = NFFT;
    holdOutput(sourcePort).df       = df;
elseif config == "line"
    holdOutput_line(sourcePort).time     = tVecs;
    holdOutput_line(sourcePort).freq     = fVecs;
    holdOutput_line(sourcePort).Z        = Z;
    holdOutput_line(sourcePort).dt       = dt;
    holdOutput_line(sourcePort).Nt       = nMax;
    holdOutput_line(sourcePort).t        = dataOut(:,1);
    holdOutput_line(sourcePort).f        = f;
    holdOutput_line(sourcePort).Nf       = NFFT;
    holdOutput_line(sourcePort).df       = df;
    holdOutput_line(sourcePort).scale    = scalingV_coeff;
end
% cd ..\..\
% save 180hybrid_wkspc.mat
% cd program\post_processing\
clear sum %where is sum redefined???????? and why???


