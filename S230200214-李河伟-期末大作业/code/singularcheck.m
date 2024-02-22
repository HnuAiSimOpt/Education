function singularcheck(kk)
% check the input matrix is singular or not
row = zeros(0);
column = zeros(0);
for count = 1:length(kk(1,:))
    maxcolumn = max(kk(:,count));
    mincolumn = min(kk(:,count));
    if maxcolumn==mincolumn
        column = [column count];
    end
    maxrow = max(kk(count,:));
    minrow = min(kk(count,:));
    if maxrow==minrow
        row = [row count];
    end
end
lenr = length(row);
lenc = length(column);
if lenr||lenc
    fprintf('\nError detected! ');
    if lenr
        fprintf('\nROW: ');
        for count = 1:lenr
            fprintf(' %d ',row(count));
        end
    end
    if lenc
        fprintf('\nCOLUMN: ');
        for count = 1:lenc
            fprintf(' %d ',column(count));
        end
        %pause(3);
    end
    fprintf('\n');
else
    fprintf('\nPASSED!\n');
end
return