function[ ] = writeBasisToFile(Z,n,c,a)
atom = '';
if Z==2
    atom = 'He';
elseif Z==4
    atom = 'Be';
elseif Z==10
    atom = 'Ne';
end

fileName = strcat('../HartreeFock/HartreeFockBases/',atom,'-STO-',num2str(n),'G');
outFile = fopen(fileName,'w'); 

fprintf(outFile,'%d %d %d %d %d\n',[Z Z/2 Z/2 Z/2 1]);
fprintf(outFile,'2 0.0 0.0 0.0\n');
fprintf(outFile,'%d\n',n);
fprintf(outFile,'0.0 0.0 0.0\n');
for i=1:n
    fprintf(outFile,'0 0 0 %15.10g %15.10g\n',[a(i) c(i)]);
end
for i=1:Z
    fprintf(outFile,'1.0\n');
end
