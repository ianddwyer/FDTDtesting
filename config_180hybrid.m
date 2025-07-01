%% CONFIGURATION FILE FOR FDTD
% 180-DEGREE HYBRID
%this file can be independently run to test PEC trace design
%the objects are useful for referencing design parameters
try
    close(5)
catch
    %
end
%use array to selecto which traces to plot
clc

low=1; high=2; 
show_trace = [1:3,4:8];%select or omit traces by listing their trace numbers here
show_substrate = true;
plotPhysical = true;

%% TUNE PARAMS, CONVERT MIL TO SAMPLE & SCALE

%discrete steps chosen to be smaller than 20/lambda and they are the size
%of the lowest common factor for the sim, which is 7mil
%unfortunately, this produces 7 mil thick PEC
cond_thk = 1.37;
LCF = 9; %lowest common factor
cm2in = 2.54; in2mil=1000; m2cm = 100;
dx = LCF/m2cm/in2mil*cm2in;
dy = dx;
dz = dx/(LCF/cond_thk); %note that discretization is 7mil but PEC is 1mil, so need to scale dz...
dl = [dx,dx,dx];
scale = round(dl(1)*m2cm*in2mil/cm2in); %scales traces, mindful of lowest common denominator


PMLdepth = 10;
airgap = 2;
noportbuffer = PMLdepth+airgap+1; %cell space from wall for PML and free space

 
%values in mils
subheight   =  50;
portlength  = 500/2;
ringlength  = qwav*2;
ringheight  = qwav;
rtracewith  =  18;
portwidth   =  45;
LRportgap   = qwav/2;

%scaling
subheight   = round( subheight  /scale )-1;
portlength  = round( portlength /scale )-1;
ringlength  = round( ringlength /scale )-1;
ringheight  = round( ringheight /scale )-1;
rtracewith  = round( rtracewith /scale )-1;
portwidth   = round( portwidth  /scale )-1;
LRportgap   = round( LRportgap  /scale )-1;

%hold for final object
linespecs(1).subheight   = subheight;
linespecs(1).portlength  = portlength;
linespecs(1).ringlength  = ringlength;
linespecs(1).ringheight  = ringheight;
linespecs(1).rtracewith  = rtracewith;
linespecs(1).portwidth   = portwidth;
linespecs(1).LRportgap   = LRportgap;

%% TEST DOMAIN CONSTRUCTION
mesh = struct([]);

%testing boundaries
mesh(1).xmax = portlength*2+ringlength;
mesh(1).ymax = portlength+ringheight+noportbuffer;
mesh(1).zmax = subheight+PMLdepth+airgap;
mesh(1).PMLdepth = PMLdepth;
mesh(1).airgap = airgap;
mesh(1).dx = dl(1);
mesh(1).dy = dl(2);
mesh(1).dz = dl(3);
mesh(1).opFreq = 3e9;
mesh(1).simtime = simtime;
mesh(1).CFLN = 0.99;

%% SUBSTRATE CONSTRUCTION
substrate = struct([]);

substrate(1).i = [1 mesh(1).xmax];
substrate(1).j = [1 mesh(1).ymax];
substrate(1).k = [1 subheight]; %actually 1:4 but costly to plot
substrate(1).epsr = 10.6;
substrate(1).sig = 0;


%% TRACE CONSTRUCTION
trace = struct([]);
%reference position is bottom left of each block
%verticle vs horizontal length is necessary to discern

%horizontal traces
n_horiz = 4; %number of horizontal traces

%lower ring bar
trace(1).jref = noportbuffer;
trace(1).iref = portlength+1;
trace(1).l = ringlength;
trace(1).w = rtracewith;
trace(1).h = subheight;

%left port
trace(2).jref = round(noportbuffer+ringheight/2)-(rtracewith+1);
trace(2).iref = 1;
trace(2).l = portlength-1;
trace(2).w = portwidth;
trace(2).h = subheight;

%right port
trace(3).jref = round(noportbuffer+ringheight/2)-(rtracewith+1);
trace(3).iref = mesh(1).xmax-portlength;
trace(3).l = portlength;
trace(3).w = portwidth;
trace(3).h = subheight;

%upper ring bar
trace(4).jref = noportbuffer+ringheight-1;
trace(4).iref = portlength;
trace(4).l = ringlength+1;
trace(4).w = rtracewith;
trace(4).h = subheight;

%verticle traces 
n_vert = 4; %number of verticle traces

%left ring bar
trace(5).jref = noportbuffer;
trace(5).iref = portlength;
trace(5).l = ringheight;
trace(5).w = rtracewith;
trace(5).h = subheight;

%right ring bar
trace(6).jref = noportbuffer;
trace(6).iref = mesh(1).xmax-portlength-rtracewith+1;
trace(6).l = ringheight;
trace(6).w = rtracewith;
trace(6).h = subheight;

%upper left port
trace(7).jref = noportbuffer+ringheight;
trace(7).iref = portlength+LRportgap-1;
trace(7).l = portlength;
trace(7).w = portwidth;
trace(7).h = subheight;

%upper right port
trace(8).jref = noportbuffer+ringheight;
trace(8).iref = mesh(1).xmax-portlength-LRportgap-rtracewith*3+1;
trace(8).l = portlength;
trace(8).w = portwidth;
trace(8).h = subheight;


%setup PEC structure from traces, this should not be modified in config
PEC = struct([]);
%horizontal sections
for num = 1:n_horiz   
    PEC(num).i = trace(num).iref+[0,trace(num).l];  PEC(num).j = trace(num).jref+[0,trace(num).w];  PEC(num).k = [trace(num).h,trace(num).h];  %left ring bar
end
%verticle sections
for num = n_horiz+1:n_horiz+n_vert
    PEC(num).i = trace(num).iref+[0,trace(num).w];  PEC(num).j = trace(num).jref+[0,trace(num).l];  PEC(num).k = [trace(num).h,trace(num).h];
end

%% PORT CONSTRUCTION
nports = 4;
NPML = 10;
%setup port blocks
port = struct([]);
%top left port (verticle)
port_marker_len = 4;
port(1).i = PEC(7).i;    port(1).j = [mesh(1).ymax-(NPML+2+port_marker_len),mesh(1).ymax-(NPML+2)];   port(1).k = substrate(1).k;
%left port (horizontal)
port(2).i = [NPML+2,NPML+port_marker_len+2];   port(2).j = PEC(2).j;   port(2).k = substrate(1).k;
%top right port (verticle)
port(3).i = PEC(8).i;    port(3).j = [mesh(1).ymax-(NPML+2+port_marker_len),mesh(1).ymax-(NPML+2)];   port(3).k = substrate(1).k;
%right port (horizontal)
port(4).i = [mesh(1).xmax-(NPML+port_marker_len+2),mesh(1).xmax-(NPML+2)];   port(4).j = PEC(3).j;   port(4).k = substrate(1).k;

%% SOURCE CONSTRUCTION
source = struct([]);

%compensate for verticle vs horizontal port directions
if sourcePort == 2 || sourcePort == 4
    source(1).i = round([mean(port(sourcePort).i),mean(port(sourcePort).i)]);
    source(1).j = port(sourcePort).j;
elseif sourcePort == 1 || sourcePort == 3
    source(1).i = port(sourcePort).i;
    source(1).j = round([mean(port(sourcePort).j),mean(port(sourcePort).j)]);
end

source(1).port = sourcePort;
source(1).k = [1,subheight];
source(1).pol = "z";
source(1).type = "SoftVoltage";
source(1).sig = "Gaussian";
source(1).tw = 1.591e-11;
source(1).nto = 5;
source(1).freq = 0;
source(1).thetInc = 0;
source(1).phiInc = 0;
source(1).resInc = 0;
source(1).amp = 1;


%% PROBE CONSTRUCTION
probe = struct([]);

shift_probe_from_source = 4;
%port 1 voltage probe
probe(1).type = "Voltage";
probe(1).axis = 'z';
probe(1).plot = 'noplot';
probe(1).i = [round(mean(port(1).i)),round(mean(port(1).i))];%[round((port(1).i(1)+port(1).i(2))/2),round((port(1).i(1)+port(1).i(2))/2)];
probe(1).j = [round(mean(port(1).j)),round(mean(port(1).j))]-shift_probe_from_source;
probe(1).k = [1,subheight];
probe(1).port = 1;

% port 2 voltage probe
probe(2).type = "Voltage";
probe(2).axis = 'z';
probe(2).plot = 'noplot';
probe(2).i = [round(mean(port(2).i)),round(mean(port(2).i))]+shift_probe_from_source;
probe(2).j = [round(mean(port(2).j)),round(mean(port(2).j))];%[round(mean(port(2).j)),round(mean(port(2).j))];
probe(2).k = [1,subheight];
probe(2).port = 2;

%port 3 voltage probe
probe(3).type = "Voltage";
probe(3).axis = 'z';
probe(3).plot = 'noplot';
probe(3).i = [round(mean(port(3).i)),round(mean(port(3).i))];%[round((port(3).i(1)+port(3).i(2))/2),round((port(3).i(1)+port(3).i(2))/2)];
probe(3).j = [round(mean(port(3).j)),round(mean(port(3).j))]-shift_probe_from_source;
probe(3).k = [1,subheight];
probe(3).port = 3;

%port 4 voltage probe
probe(4).type = "Voltage";
probe(4).axis = 'z';
probe(4).plot = 'noplot';
probe(4).i = [round(mean(port(4).i)),round(mean(port(4).i))]-shift_probe_from_source;
probe(4).j = [round(mean(port(4).j)),round(mean(port(4).j))];%[round((port(4).j(1)+port(4).j(2))/2),round((port(4).j(1)+port(4).j(2))/2)];
probe(4).k = [1,subheight];
probe(4).port = 4;



%Current Probe 1 for Impedance
probe(13).type = "Current";
probe(13).plot = 'noplot';
probe(13).k = [subheight,subheight];
probe(13).port = sourcePort;

%Current Probe 2 for Impedance
probe(14).type = "Current";
probe(14).plot = 'noplot';
probe(14).k = [subheight,subheight];
probe(14).port = sourcePort;

%Current Probe 3 for Impedance
probe(5).type = "Current";
probe(5).plot = 'noplot';
probe(5).axis = 'y';
probe(5).i = [round(mean(port(1).i))-1,round(mean(port(1).i))-1];
probe(5).j = [round(mean(port(1).j)),round(mean(port(1).j))+1]-shift_probe_from_source;
probe(5).k = [subheight,subheight];
probe(5).port = 1;

%Current Probe 4 for Impedance
probe(6).type = "Current";
probe(6).plot = 'noplot';
probe(6).axis = 'y';
probe(6).i = [round(mean(port(1).i))+1,round(mean(port(1).i))+1];
probe(6).j = [round(mean(port(1).j)),round(mean(port(1).j))+1]-shift_probe_from_source;
probe(6).k = [subheight,subheight];
probe(6).port = 1;

%Current Probe 5 for Impedance
probe(7).type = "Current";
probe(7).plot = 'noplot';
probe(7).axis = 'x';
probe(7).i = [round(mean(port(2).i)),round(mean(port(2).i))+1]+shift_probe_from_source;
probe(7).j = [round(mean(port(2).j))-1,round(mean(port(2).j))-1];
probe(7).k = [subheight,subheight];
probe(7).port = 2;

%Current Probe 6 for Impedance
probe(8).type = "Current";
probe(8).plot = 'noplot';
probe(8).axis = 'x';
probe(8).i = [round(mean(port(2).i)),round(mean(port(2).i))+1]+shift_probe_from_source;
probe(8).j = [round(mean(port(2).j))+1,round(mean(port(2).j))+1];
probe(8).k = [subheight,subheight];
probe(8).port = 2;

%Current Probe 7 for Impedance
probe(9).type = "Current";
probe(9).plot = 'noplot';
probe(9).axis = 'y';
probe(9).i = [round(mean(port(3).i))-1,round(mean(port(3).i))-1];
probe(9).j = [round(mean(port(3).j)),round(mean(port(3).j))+1]-shift_probe_from_source;
probe(9).k = [subheight,subheight];
probe(9).port = 3;

%Current Probe 8 for Impedance
probe(10).type = "Current";
probe(10).plot = 'noplot';
probe(10).axis = 'y';
probe(10).i = [round(mean(port(3).i))+1,round(mean(port(3).i))+1];
probe(10).j = [round(mean(port(3).j)),round(mean(port(3).j))+1]-shift_probe_from_source;
probe(10).k = [subheight,subheight];
probe(10).port = 3;


%Current Probe 7 for Impedance
probe(11).type = "Current";
probe(11).plot = 'noplot';
probe(11).axis = 'x';
probe(11).i = [round(mean(port(4).i)),round(mean(port(4).i))+1]-shift_probe_from_source;
probe(11).j = [round(mean(port(4).j))-1,round(mean(port(4).j))-1];
probe(11).k = [subheight,subheight];
probe(11).port = 4;

%Current Probe 8 for Impedance
probe(12).type = "Current";
probe(12).plot = 'noplot';
probe(12).axis = 'x';
probe(12).i = [round(mean(port(4).i)),round(mean(port(4).i))+1]-shift_probe_from_source;
probe(12).j = [round(mean(port(4).j))+1,round(mean(port(4).j))+1];
probe(12).k = [subheight,subheight];
probe(12).port = 4;

%NOTE THAT SOURCE CURRENT PROBES FOR Z NEED TO BE LAST 2 OF STRUCT TO WORK WITH
%CODE

%compensate for verticle vs horizontal port directions
shift = 0; %fix for offset current probe
if sourcePort == 2 || sourcePort == 4
    probe(13).axis = 'x';
    probe(14).axis = 'x';
    if sourcePort==2
        shift=0;
        probe(13).i = [round(mean(source(1).i)),round(mean(source(1).i))+1]+shift+shift_probe_from_source;
        probe(13).j = [round(mean(source(1).j))-1,round(mean(source(1).j))-1];
        probe(14).i = [round(mean(source(1).i)),round(mean(source(1).i))+1]+shift+shift_probe_from_source;
        probe(14).j = [round(mean(source(1).j))+1,round(mean(source(1).j))+1];
    end
    if sourcePort==4
        probe(13).i = [round(mean(source(1).i)),round(mean(source(1).i))+1]+shift-shift_probe_from_source;
        probe(13).j = [round(mean(source(1).j))-1,round(mean(source(1).j))-1];
        probe(14).i = [round(mean(source(1).i)),round(mean(source(1).i))+1]+shift-shift_probe_from_source;
        probe(14).j = [round(mean(source(1).j))+1,round(mean(source(1).j))+1];
    end
elseif sourcePort == 1 || sourcePort == 3
    probe(13).axis = 'y';
    probe(14).axis = 'y';
    probe(13).i = [round(mean(source(1).i))-1,round(mean(source(1).i))-1];
    probe(13).j = [round(mean(source(1).j)),round(mean(source(1).j))+1]-shift_probe_from_source;
    probe(14).i = [round(mean(source(1).i))+1,round(mean(source(1).i))+1];
    probe(14).j = [round(mean(source(1).j)),round(mean(source(1).j))+1]-shift_probe_from_source;
end

%% PHYSICAL DOMAIN VISUALIZATION

sub_block = zeros([3 mesh(1).ymax*mesh(1).zmax*mesh(1).xmax]);
port_block = zeros([3 mesh(1).ymax*mesh(1).zmax*mesh(1).xmax]);
pec_block = zeros([3 mesh(1).ymax*mesh(1).zmax*mesh(1).xmax]);
count = 1;
if plotPhysical
    if show_substrate
        for num = 1:1
            for kk = substrate(num).k(low):substrate(num).k(high)
                for jj = substrate(num).j(low):substrate(num).j(high)
                    for ii = substrate(num).i(low):substrate(num).i(high)
                        sub_block(1,count) = ii; sub_block(2,count) = jj; sub_block(3,count) = kk; 
                        count=count+1;
                    end
                end
            end
        end
    end
    count = 1;
    for num = 1:4 %total ports
        for kk = port(num).k(low):port(num).k(high)
            for jj = port(num).j(low):port(num).j(high)
                for ii = port(num).i(low):port(num).i(high)
                    port_block(1,count) = ii; port_block(2,count) = jj; port_block(3,count) = kk; 
                    count=count+1;
                end
            end
        end
    end
    count = 1;
    for seg = show_trace
        for kk = PEC(seg).k(low):PEC(seg).k(high)
            for jj = PEC(seg).j(low):PEC(seg).j(high)
                for ii = PEC(seg).i(low):PEC(seg).i(high)
                    pec_block(1,count) = ii; pec_block(2,count) = jj; pec_block(3,count) = kk; 
                    count=count+1;
                end
            end
        end
    end

    figure('name','PCB')
    scatter3(sub_block(1,:),sub_block(2,:),sub_block(3,:),'g','o');
    hold on
    scatter3(pec_block(1,:),pec_block(2,:),pec_block(3,:),'b','.');
    scatter3(port_block(1,:),port_block(2,:),port_block(3,:),'r','o');
    ylabel('j(index)')
    xlabel('i(index)')
    zlabel('k(index)')
    view([0 90])
    title("PHYSICAL CONSTRUCTION")
    legend("substrate","PEC","port")
    
    %if x>y, use to center board on display
    spacexy =(mesh(1).xmax-mesh(1).ymax)/2;
    %if x>z, use to center board on display
    spacexz = (mesh(1).xmax-mesh(1).zmax)/2;

    ylim([-spacexy mesh(1).ymax+spacexy])
    xlim([-10 mesh(1).xmax+spacexy])
    zlim([-spacexz mesh(1).zmax+spacexz]) %increasing z will cause it to appear more dense
    scale_window = (mesh(1).ymax+2*spacexy)/(mesh(1).xmax+2*spacexy);
    grid on
    set(gcf,'Position',[100 100 round(900) round(900*scale_window)])
    %label locations 
    text(probe(7).i,probe(7).j,probe(7).k,"I-probe1","FontWeight",'bold',FontSize=8)
    text(probe(8).i,probe(8).j,probe(8).k,"I-probe2","FontWeight",'bold',FontSize=8)
    text(probe(9).i,probe(9).j,probe(9).k,"I-probe1","FontWeight",'bold',FontSize=8)
    text(probe(10).i,probe(10).j,probe(10).k,"I-probe2","FontWeight",'bold',FontSize=8)
    text(probe(11).i,probe(11).j,probe(11).k,"I-probe1","FontWeight",'bold',FontSize=8)
    text(probe(12).i,probe(12).j,probe(12).k,"I-probe2","FontWeight",'bold',FontSize=8)
    text(probe(13).i,probe(13).j,probe(13).k,"I-probe1","FontWeight",'bold',FontSize=8)
    text(probe(14).i,probe(14).j,probe(14).k,"I-probe2","FontWeight",'bold',FontSize=8)
    text(probe(6).i,probe(6).j,probe(6).k,"I-probe1","FontWeight",'bold',FontSize=8)
    text(probe(3).i,probe(3).j,probe(3).k,"port(3)"+newline+"V-probe(2)","FontWeight",'bold',FontSize=8)
    text(probe(5).i,probe(5).j,probe(5).k,"I-probe2","FontWeight",'bold',FontSize=8)
    text(probe(4).i,probe(4).j,probe(4).k,"port(4)"+newline+"V-probe(4)","FontWeight",'bold',FontSize=8)
    text(probe(2).i,probe(2).j,probe(2).k,"port(2)"+newline+"V-probe(5)","FontWeight",'bold',FontSize=8)
    text(probe(1).i,probe(1).j,probe(1).k,"port(1)"+newline+"V-probe(6)","FontWeight",'bold',FontSize=8)
    text(source(1).i,source(1).j,source(1).k,"source","FontWeight",'bold',FontSize=8)
end
clear sub_block pec_block port_block

%% fixes mimsuse of max earlier in script
%lengths were extended by 1 to manage overlaps to ensure intended distance
mesh(1).xmax = mesh(1).xmax;
mesh(1).ymax = mesh(1).ymax;
mesh(1).zmax = mesh(1).zmax;
