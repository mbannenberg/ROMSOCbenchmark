% ------------------------------------------------------------------------------
% Netlist writing function
%
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

fileID = fopen('Circuits/multirate_example_long_nonlin.cir','w');
fprintf(fileID,'V%.d %s %d sin 5 40\n',1,'SUB',1);
for i = 1:2
    fprintf(fileID,'D%.d %d %d 1e-12\n',i,i,i+1);
    fprintf(fileID,'R%.d %d %s 1e3\n',i,i+1,'SUB');
    fprintf(fileID,'C%.d %d %s 10e-6\n',i,i+1,'SUB');
end
for i = 3:3000
    fprintf(fileID,'R%.d %d %d %de4\n',i,3,i+1,i);
    fprintf(fileID,'D%.d %d %d 1e-12\n',i,i+1,i+2);
    fprintf(fileID,'C%.d %d %s 10e-6\n',i,i+2,'SUB');
end 
    
fclose(fileID);