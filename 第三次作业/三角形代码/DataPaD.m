function [x,y,z,M,tx,ty,tz,Mt] = DataPaD(Disp,nns,ndpn)

x = [];
y = [];
z = [];
tx = [];
ty = [];
tz = [];
for count = 1:nns
     x(count) = Disp(ndpn*(count-1)+1);
     y(count) = Disp(ndpn*(count-1)+2);
     if ndpn>2
        z(count) = Disp(ndpn*(count-1)+3);
     else
         z(count)=0;
     end
     if ndpn>3
         tx(count) = Disp(ndpn*(count-1)+4);
         ty(count) = Disp(ndpn*(count-1)+5);
         tz(count) = Disp(ndpn*(count-1)+6);
     else
         tx(count) = 0;
         ty(count) = 0;
         tz(count) = 0;
     end
end
all = {x,y,z,tx,ty,tz};
for count=1:length(all)
    if max(all{count})==min(all{count})
        fprintf('Disp Err:%d\n',count);
    end
end
for count = 1:length(z)
     M(count) = sqrt( y(count)^2+ x(count)^2+ z(count)^2);
    Mt(count) = sqrt(ty(count)^2+tx(count)^2+tz(count)^2);
end




