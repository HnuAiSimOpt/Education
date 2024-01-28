function [node,element]=Rect_mesh(a,b,Nx,Ny,etype)
% To acquire mesh for the rectangle plate
% a,b---length of two sides
% Nx,Ny---number of division in each direction
% etype---element type, 'tri' and 'quad'
% 矩形板的网格生成函数，a,b,为矩形板的长宽，Nx, Ny 为矩形板在x,y方向划分的单元数
% etype 为单元类型，取'tri'是为三角形单元，取'quad'时为四边形单元


x=linspace(0,a,Nx+1);
y=linspace(0,b,Ny+1);

node=zeros((Nx+1)*(Ny+1),2);
count=0;   % count for nodes
for ix=1:Nx+1
    for iy=1:Ny+1
        count=count+1;
        node(count,:)=[x(ix),y(iy)];
    end
end

% define a function for nodes
f=@(set) (set(1)-1)*(Ny+1)+set(2);

if strcmp(etype,'tri')
    element=zeros(2*Nx*Ny,3);
    count=0; % count for elements
    for ix=1:Nx
        for iy=1:Ny
            count=count+1;
            element(count,:)=[f([ix,iy]),f([ix+1,iy]),f([ix+1,iy+1])];
            count=count+1;
            element(count,:)=[f([ix,iy]),f([ix+1,iy+1]),f([ix,iy+1])];
        end
    end
elseif strcmp(etype,'quad')
    element=zeros(Nx*Ny,4);
    count=0; % count for elements
    for ix=1:Nx
        for iy=1:Ny
            count=count+1;
            element(count,:)=[f([ix,iy]),f([ix+1,iy]),f([ix+1,iy+1]),f([ix,iy+1])];
        end
    end
end