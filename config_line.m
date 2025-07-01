%% CONFIGURATION FILE FOR FDTD
% microstrip transmission line

%this file can be independently run to test PEC trace design
%the objects are useful for referencing design parameters

%use array to selecto which traces to plot
clc
low=1; high=2; 
show_trace = 1:1;
show_substrate = true;
plotPhysical = true;


%% mil conversion to sample

%matches 180 hybrid for ref
cond_thk = 1.39;
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

%% TEST DOMAIN CONSTRUCTION
mesh = struct([]);

%testing boundaries
mesh(1).xmax = round((portlength*2+ringlength));
mesh(1).ymax = round((portlength+ringheight+noportbuffer));
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
substrate(1).k = [1 subheight];
substrate(1).epsr = 10.6;
substrate(1).sig = 0;


%% TRACE CONSTRUCTION
trace = struct([]);
%reference position is bottom left of each block
%verticle vs horizontal length is necessary to discern

%horizontal traces
n_horiz = 0; %number of horizontal PEC traces

%lower ring bar
% trace(1).jref = noportbuffer;
% trace(1).iref = linelength;
% trace(1).l = ringlength;
% trace(1).w = ring_trace_width;
% trace(1).h = sub_height;

%verticle traces 
n_vert = 1; %number of verticle traces

%upper left port
trace(1).jref = 1;
trace(1).iref = round(mesh(1).xmax/2-portwidth/2);
trace(1).l = mesh(1).ymax-1;
trace(1).w = portwidth;
trace(1).h = subheight;


%setup PEC structure from traces, this should not be modified in config
PEC = struct([]);

%verticle sections
for num = n_horiz+1:n_horiz+n_vert
    PEC(num).i = trace(num).iref+[0,trace(num).w];  PEC(num).j = trace(num).jref+[0,trace(num).l];  PEC(num).k = [trace(num).h,trace(num).h];
end

%% PORT CONSTRUCTION
nports = 2;
PMLdepth = 10;
airgap = 3;
noportbuffer = 2*airgap+PMLdepth;
%setup port blocks
port = struct([]);
%top left port (verticle)
port_depth = 4;
port(1).i = PEC(1).i;    port(1).j = [noportbuffer,noportbuffer+port_depth];   port(1).k = substrate(1).k;
port(2).i = PEC(1).i;    port(2).j = [trace(1).l-noportbuffer-port_depth,trace(1).l-noportbuffer];   port(2).k = substrate(1).k;


%% SOURCE CONSTRUCTION
source = struct([]);


source(1).i = port(sourcePort).i;
source(1).j = round([mean(port(sourcePort).j),mean(port(sourcePort).j)]);

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
if sourcePort==2
    shift_probe_from_source = -shift_probe_from_source ;
end

%voltage probes
probe(1).type = "Voltage";
probe(1).axis = 'z';
probe(1).plot = 'noplot';
probe(1).i = [round(mean(port(sourcePort).i)),round(mean(port(sourcePort).i))];%[round((port(sourcePort).i(1)+port(sourcePort).i(2))/2),round((port(sourcePort).i(1)+port(sourcePort).i(2))/2)];
probe(1).j = [round(mean(port(sourcePort).j)),round(mean(port(sourcePort).j))]+shift_probe_from_source;
probe(1).k = [1,subheight];
probe(1).port = sourcePort;

probe(2).type = "Voltage";
probe(2).axis = 'z';
probe(2).plot = 'noplot';
probe(2).i = [round(mean(port(sourcePort).i)),round(mean(port(sourcePort).i))];%[round((port(sourcePort).i(1)+port(sourcePort).i(2))/2),round((port(sourcePort).i(1)+port(sourcePort).i(2))/2)];
probe(2).j = [round(mean(port(sourcePort).j)),round(mean(port(sourcePort).j))];
probe(2).k = [1,subheight];
probe(2).port = sourcePort;

%NOTE THAT CURRENT PROBES FOR Z NEED TO BE LAST 2 OF STRUCT TO WORK WITH
%CODE
%current probes 
probe(3).type = "Current";
probe(3).axis = 'y';
probe(3).plot = 'noplot';
probe(3).i = [round(mean(port(sourcePort).i))-1,round(mean(port(sourcePort).i))-1];
probe(3).j = [round(mean(port(sourcePort).j)),round(mean(port(sourcePort).j))+1];
probe(3).k = [subheight,subheight];
probe(3).port = sourcePort;
%current probes
probe(4).type = "Current";
probe(4).axis = 'y';
probe(4).plot = 'noplot';
probe(4).i = [round(mean(port(sourcePort).i))+1,round(mean(port(sourcePort).i))+1];
probe(4).j = [round(mean(port(sourcePort).j)),round(mean(port(sourcePort).j))+1];
probe(4).k = [subheight,subheight];
probe(4).port = sourcePort;

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
    max(sub_block,[],[ 2 3])
    count = 1;
    for num = 1:length(port) %total ports
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
    title("PHYSICAL CONSTRUCTION")
    legend("substrate","PEC","port")
    
    ylim([0 mesh(1).xmax+10])
    xlim([-10 mesh(1).xmax+10])
    zlim([-10 mesh(1).xmax+10]) %increasing z will cause it to appear more dense

    grid on
    set(gcf,'Position',[100 100 round(900*3/4) round(720*3/4)])
    %label locations 
    text(probe(length(probe)-1).i,probe(length(probe)-1).j,probe(length(probe)-1).k,"I-probe1","FontWeight",'bold',FontSize=8)
    text(probe(1).i,probe(1).j,probe(1).k,"port(1)"+newline+"V-probe(1)","FontWeight",'bold',FontSize=8)
    text(probe(length(probe)).i,probe(length(probe)).j,probe(length(probe)).k,"I-probe2","FontWeight",'bold',FontSize=8)
    text(source(1).i,source(1).j,source(1).k,"source","FontWeight",'bold',FontSize=8)
end
clear sub_block pec_block port_block

%% STRANGE INDEXING FIX

mesh(1).xmax = mesh(1).xmax;
mesh(1).ymax = mesh(1).ymax;
mesh(1).zmax = mesh(1).zmax;
