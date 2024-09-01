clear;
clc;

nr = 2;
nctrl = 3;

i = 0;
for r1 = 1:nr
    for c1 = 1:nctrl
        i = i+1;
        j = 0;
        
        for r = 1:nr
            for c = 1:nctrl
                j = j+1;
                disp([num2str(i) num2str(j) ':' num2str(r) num2str(r1) num2str(c) num2str(c1)])
                %matrix(i,j) = num2str(r)+num2str(r1)+num2str(c)+num2str(c1);
                %disp(num2str(r)+num2str(r1)+num2str(c)+num2str(c1))
            end
        end
    end
end

disp('new')
%disp(matrix)
i = 0;
for r = 1:nr
    for c1 = 1:nctrl
    i =  i + 1;
    j = 0;
        for r1 = 1:nr
            for c = 1:nctrl
                j = j + 1;
                disp([num2str(i) num2str(j) ':' num2str(r) num2str(r1) num2str(c) num2str(c1)])
            end
        end
    end
end


