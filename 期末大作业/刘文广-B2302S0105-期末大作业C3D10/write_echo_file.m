function linestr=write_echo_file(fid_inp)
c_Tab=char(9);
readOK=0;

while readOK==0    
    linestr = fgets(fid_inp);    
    lengstr=length(linestr);
    for i=1:lengstr
        c=linestr(i:i);
        if (c== '%' | c=='#'| c=='*'| c=='$')
            break
        elseif (c==' ' | c==c_Tab) 
        else
            readOK=1;
        end
    end
end