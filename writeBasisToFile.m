function[ ] = writeBasisToFile(Z,n,c,a,basisSize)
atom = '';
if Z==2
    atom = 'He';
elseif Z==4
    atom = 'Be';
elseif Z==10
    atom = 'Ne';
else 
    exit(1);
end

fileName = strcat('../HartreeFock/HartreeFockBases/',atom,'-STO-',num2str(n),'G');
outFile = fopen(fileName,'w');

fprintf(outFile,'%d %d %d %d %d\n',[Z Z/2 Z/2 Z/2 1]);
fprintf(outFile,'%d 0.0 0.0 0.0\n',Z);

fprintf(outFile,'%d\n',n);
fprintf(outFile,'0.0 0.0 0.0\n');
if Z==2
    for i=1:n
        fprintf(outFile,'0 0 0 %15.10g %15.10g\n',[a(i) c(i)]);
    end
else
    if Z>3
        b=1;  %   1s
        for i=1:n
            fprintf(outFile,'0 0 0 %15.10g %15.10g\n',[a(b,i) c(b,i)]);
        end
        
        fprintf(outFile,'%d\n',n);
        fprintf(outFile,'0.0 0.0 0.0\n');
        b=2;  %   2s
        for i=1:n
            fprintf(outFile,'0 0 0 %15.10g %15.10g\n',[a(b,i) c(b,i)]);
        end
    end
    if Z>5
        fprintf(outFile,'%d\n',n);
        fprintf(outFile,'0.0 0.0 0.0\n');
        b=3;  %   2px
        for i=1:n
            fprintf(outFile,'1 0 0 %15.10g %15.10g\n',[a(b,i) c(b,i)]);
        end
        fprintf(outFile,'%d\n',n);
        fprintf(outFile,'0.0 0.0 0.0\n');
        b=3;  %   2py
        for i=1:n
            fprintf(outFile,'0 1 0 %15.10g %15.10g\n',[a(b,i) c(b,i)]);
        end
        fprintf(outFile,'%d\n',n);
        fprintf(outFile,'0.0 0.0 0.0\n');
        b=3;  %   2pz
        for i=1:n
            fprintf(outFile,'0 0 1 %15.10g %15.10g\n',[a(b,i) c(b,i)]);
        end
    end
end

if Z==2
    fprintf(outFile,'1.0\n');
    fprintf(outFile,'1.0\n');
elseif Z==4
    fprintf(outFile,'1.0 0.0\n');
    fprintf(outFile,'0.0 1.0\n');
    fprintf(outFile,'1.0 0.0\n');
    fprintf(outFile,'0.0 1.0\n');
else
    fprintf(outFile,'1.0 0.0 0.0 0.0 0.0\n');
    fprintf(outFile,'0.0 1.0 0.0 0.0 0.0\n');
    fprintf(outFile,'0.0 0.0 1.0 0.0 0.0\n');
    fprintf(outFile,'0.0 0.0 0.0 1.0 0.0\n');
    fprintf(outFile,'0.0 0.0 0.0 0.0 1.0\n');
    fprintf(outFile,'1.0 0.0 0.0 0.0 0.0\n');
    fprintf(outFile,'0.0 1.0 0.0 0.0 0.0\n');
    fprintf(outFile,'0.0 0.0 1.0 0.0 0.0\n');
    fprintf(outFile,'0.0 0.0 0.0 1.0 0.0\n');
    fprintf(outFile,'0.0 0.0 0.0 0.0 1.0\n');
end
end
