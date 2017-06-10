data = importdata('OCSEt_140000_mat.txt');

Xcoords = data(:,2);
Ycoords = data(:,3);
Zcoords = data(:,4);

% Xcoords = data(:,1);
% Ycoords = data(:,2);
% Zcoords = data(:,3);

%TRI = delaunay(Xcoords,Ycoords,Zcoords);

%size(TRI);

%scatter3(Xcoords,Ycoords,Zcoords);
%trimesh(TRI,Xcoords,Ycoords,Zcoords);

DT = delaunayTriangulation(Xcoords,Ycoords,Zcoords)
pts = DT.Points;
xpts = pts(:,1);
ypts = pts(:,2);
zpts = pts(:,3);

%scatter3(xpts,ypts,zpts)
%tetramesh(DT)

con = DT.ConnectivityList;

con(1,:)
[m,n] = size(con)


for i = 1:m
    for j = 1:n
        
        
    end   
end
