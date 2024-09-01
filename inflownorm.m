load VBZ_PSI.dat;
rad=8.178;
asr=3*pi/180;
m2f=3.28084;

norm=0;

k=0;
for p=1:36
    for s=1:32
        k=k+1;
        inflownw(s,p)=VBZ_PSI(k,6);
        norm=norm+inflownw(s,p);
    end
end    

normmps=-norm/(32*36)

% normfps=normmps*m2f

normnd=normmps/(124.6*1.676)