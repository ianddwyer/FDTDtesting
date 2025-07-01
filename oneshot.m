%init object structures for outputs
holdOutput = struct([]);
holdOutput_line = struct([]);
fVecs = struct([]);
tVecs = struct([]);
linespecs = struct([]); %holds line dimensions
wrkspc_name = sprintf(string(config)+"_wkspc.mat");  
cd program;
Main3DFDTD_Accelerated; %runs main with new configuration for sourcePort 
cd ..;