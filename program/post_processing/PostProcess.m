%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Post Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A lot still TBD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output & Plot designated output values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_n = dt*(0:nMax-1);
time_nph = time_n+0.5*dt;
%  Loop over the output quantities:

for iout = 1:numOutputQty
    
    if OutputType(iout) <= 8
        switch OutputType(iout)
            case 1
                titleStr = ['Ex Probe ' num2str(iout,'%d')];
                ylabelStr = 'Ex (V/m)';
                filename = ['ExProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_n';
            case 2
                titleStr = ['Ey Probe ' num2str(iout,'%d')];
                ylabelStr = 'Ey (V/m)';
                filename = ['EyProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_n';
            case 3
                titleStr = ['Ez Probe ' num2str(iout,'%d')];
                ylabelStr = 'Ez (V/m)';
                filename = ['EzProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_n';
            case 4
                titleStr = ['Hx Probe ' num2str(iout,'%d')];
                ylabelStr = 'Hx (A/m)';
                filename = ['HxProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_nph';
            case 5
                titleStr = ['Hy Probe ' num2str(iout,'%d')];
                ylabelStr = 'Hy (A/m)';
                filename = ['HyProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_nph';
            case 6
                titleStr = ['Hz Probe ' num2str(iout,'%d')];
                ylabelStr = 'Hz (A/m)';
                filename = ['HzProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_nph';
            case 7
                titleStr = ['Voltage Probe ' num2str(iout,'%d') ];
                ylabelStr = 'V (V)';
                filename = ['VProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_n';
            case 8
                titleStr = ['Current Probe ' num2str(iout,'%d') ];
                ylabelStr = 'I (A)';
                filename = ['IProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_nph';
                
        end
        
        % plot if requested
        if outputPlot(iout) > 0
            figure;
            plot(dataOut(:,1),OutputValue(:,iout),'Linewidth',2);

            title(titleStr);
            xlabel('t (s)');
            ylabel(ylabelStr);
            grid on
        end

        % write data to txt file:
        [fid, msg] = fopen(filename,'w');

        dataOut(:,2) = OutputValue(:,iout);
        if fid > 0
            for no = 1:nMax
                fprintf(fid,'%d %d \n',dataOut(no,1),dataOut(no,2));
            end
            fclose(fid);
        else
            msg
        end
    end
    
end

FFT_storeData;
% save(wrkspc_name); %save before using the loads in networks params and moving on to next source loop in case of break
