figure(6); 
plot(holdOutput_line(1).t,holdOutput_line(1).time(1).data)
xlabel("time[s]")
ylabel("voltage[V]")
title("Test Line Voltage for Reference Values")
grid on

figure(8);
for ii = 1:4
    subplot(220+ii)
    plot(holdOutput(ii).f,holdOutput(ii).Z)
    title("180hybrid Volage from Port 1 to Port "+ii)
    grid on
    xlim([0 6]*1e9)
    xlabel("time[s]")
    ylabel("voltage[V]")
end
try
    close(9)
catch
end
figure(9);
legend_names = ["","","",""];
for jj = 1:length(sourceport)
    subplot(220+jj)
    for ii = 1:length(port)
        plot(f,abs(real(holdOutput(ii).Zports(jj,:)))); xlim([0 6e9])
        hold on;        
        legend_names(ii) = "Z"+string(jj)+string(ii);    
    end
    title("Impedance Over Frequency: Source at "+string(jj))
    legend(legend_names)
    grid on
end

close(8)
figure(8);
plot(holdOutput(1).f,holdOutput_line(1).Z)
title("180hybrid Reference Line Impedence")
grid on
ylim ([0 100])
xlim([0 6]*1e9)
xlabel("frequency[Hz]")
ylabel("Impedance[\Omega]")
%% CORRECT FFT SCALING FOR VOLTAGE!!!!!
figure(12);plot(1:NFFT,20*log10(abs(fft(OutputValue(:,1),NFFT)./(NFFT*dz))));xlim([0 480*2])
