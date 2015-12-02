clc;
clear all;

load('centaur_HKS.mat');
infile='centaur1.off';

fid=fopen(infile);
fgetl(fid);
nos=fscanf(fid, '%d %d %d',[3 1]);
nopts=nos(1);
notrg=nos(2);

coord=fscanf(fid, '%g %g %g',[3 nopts]);
coord=coord';

triang=fscanf(fid, '%d %d %d %d',[4 notrg]);
triang=triang';
triang=triang(:,2:4)+1;

fclose(fid);

Edist= Euclidean_distance(coord,coord);

G=sparse(nopts,nopts);

for i=1:notrg
    p1=triang(i,1);
    p2=triang(i,2);
    p3=triang(i,3);
    G(p1,p2)= Edist(p1,p2);
    G(p2,p1)= Edist(p2,p1);
    G(p1,p3)= Edist(p1,p3);
    G(p3,p1)= Edist(p3,p1);
    G(p2,p3)= Edist(p2,p3);
    G(p3,p2)= Edist(p3,p2);
end
Gdist=graphallshortestpaths(G);

inde=1:1:nopts;
temp_hks=centaur1_hks';

temclust=mean(temp_hks);

[ma,ind]=max(temclust);

center(1,1)=inde(ind);


for b=1:1:9
for i=1:1:b
   tempdist(i,:)=Gdist(center(b,1),:);
end
farth=center;
[far,fartem]=sort(sum(tempdist),'descend');
a=1;
while(any(fartem(1,a)==farth))
    a=a+1;
end
farth(b+1,1)=fartem(1,a);
tempdist=[];

    for i=1:1:nopts
        disMin=inf;
        indexPointMin=0;
        for j=1:1:size(farth,1)
            if Gdist(i,farth(j,1))<disMin
                indexPointMin=j;
                disMin=Gdist(i,farth(j));
            end
        end
        f_vec(1,i)=indexPointMin;
    end
    ind1=find(f_vec==size(farth,1));
    temp_hks=[];
    for i=1:1:size(ind1,2)
        temp_hks(i,:)=centaur1_hks(ind1(1,i),:);
    end
    temp_hks=temp_hks';
    inde=ind1;

    temclust=mean(temp_hks);

    [ma,ind]=max(temclust);

    center(b+1,1)=inde(ind);

end

for i=1:1:nopts
        disMin=inf;
        indexPointMin=0;
        for j=1:1:size(center,1)
            if Gdist(i,center(j,1))<disMin
                indexPointMin=j;
                disMin=Gdist(i,center(j,1));
            end
        end
        f_vec(1,i)=indexPointMin;
end

f_vec1=f_vec;

fprintf('Program paused. Press enter to apply k center to centaur 2\n');
pause;
center=[];
infile='centaur2.off';

fid=fopen(infile);
fgetl(fid);
nos=fscanf(fid, '%d %d %d',[3 1]);
nopts=nos(1);
notrg=nos(2);

coord=fscanf(fid, '%g %g %g',[3 nopts]);
coord=coord';

triang=fscanf(fid, '%d %d %d %d',[4 notrg]);
triang=triang';
triang=triang(:,2:4)+1;

fclose(fid);

Edist= Euclidean_distance(coord,coord);

G=sparse(nopts,nopts);

for i=1:notrg
    p1=triang(i,1);
    p2=triang(i,2);
    p3=triang(i,3);
    G(p1,p2)= Edist(p1,p2);
    G(p2,p1)= Edist(p2,p1);
    G(p1,p3)= Edist(p1,p3);
    G(p3,p1)= Edist(p3,p1);
    G(p2,p3)= Edist(p2,p3);
    G(p3,p2)= Edist(p3,p2);
end
Gdist=graphallshortestpaths(G);

inde=1:1:nopts;
temp_hks=centaur2_hks';

temclust=mean(temp_hks);

[ma,ind]=max(temclust);

center(1,1)=inde(ind);


for b=1:1:9
for i=1:1:b
   tempdist(i,:)=Gdist(center(b,1),:);
end
farth=center;
[far,fartem]=sort(sum(tempdist),'descend');
a=1;
while(any(fartem(1,a)==farth))
    a=a+1;
end
farth(b+1,1)=fartem(1,a);
tempdist=[];

    for i=1:1:nopts
        disMin=inf;
        indexPointMin=0;
        for j=1:1:size(farth,1)
            if Gdist(i,farth(j,1))<disMin
                indexPointMin=j;
                disMin=Gdist(i,farth(j));
            end
        end
        f_vec(1,i)=indexPointMin;
    end
    ind1=find(f_vec==size(farth,1));
    temp_hks=[];
    for i=1:1:size(ind1,2)
        temp_hks(i,:)=centaur2_hks(ind1(1,i),:);
    end
    temp_hks=temp_hks';
    inde=ind1;

    temclust=mean(temp_hks);

    [ma,ind]=max(temclust);

    center(b+1,1)=inde(ind);

end

for i=1:1:nopts
        disMin=inf;
        indexPointMin=0;
        for j=1:1:size(center,1)
            if Gdist(i,center(j,1))<disMin
                indexPointMin=j;
                disMin=Gdist(i,center(j,1));
            end
        end
        f_vec(1,i)=indexPointMin;
end

f_vec2=f_vec;

fprintf('Program paused. Press enter to apply k center to the centaur 3\n');
pause;

center=[];
infile='centaur3.off';

fid=fopen(infile);
fgetl(fid);
nos=fscanf(fid, '%d %d %d',[3 1]);
nopts=nos(1);
notrg=nos(2);

coord=fscanf(fid, '%g %g %g',[3 nopts]);
coord=coord';

triang=fscanf(fid, '%d %d %d %d',[4 notrg]);
triang=triang';
triang=triang(:,2:4)+1;

fclose(fid);

Edist= Euclidean_distance(coord,coord);

G=sparse(nopts,nopts);

for i=1:notrg
    p1=triang(i,1);
    p2=triang(i,2);
    p3=triang(i,3);
    G(p1,p2)= Edist(p1,p2);
    G(p2,p1)= Edist(p2,p1);
    G(p1,p3)= Edist(p1,p3);
    G(p3,p1)= Edist(p3,p1);
    G(p2,p3)= Edist(p2,p3);
    G(p3,p2)= Edist(p3,p2);
end
Gdist=graphallshortestpaths(G);

inde=1:1:nopts;
temp_hks=centaur3_hks';

temclust=mean(temp_hks);

[ma,ind]=max(temclust);

center(1,1)=inde(ind);


for b=1:1:9
for i=1:1:b
   tempdist(i,:)=Gdist(center(b,1),:);
end
farth=center;
[far,fartem]=sort(sum(tempdist),'descend');
a=1;
while(any(fartem(1,a)==farth))
    a=a+1;
end
farth(b+1,1)=fartem(1,a);
tempdist=[];

    for i=1:1:nopts
        disMin=inf;
        indexPointMin=0;
        for j=1:1:size(farth,1)
            if Gdist(i,farth(j,1))<disMin
                indexPointMin=j;
                disMin=Gdist(i,farth(j));
            end
        end
        f_vec(1,i)=indexPointMin;
    end
    ind1=find(f_vec==size(farth,1));
    temp_hks=[];
    for i=1:1:size(ind1,2)
        temp_hks(i,:)=centaur3_hks(ind1(1,i),:);
    end
    temp_hks=temp_hks';
    inde=ind1;

    temclust=mean(temp_hks);

    [ma,ind]=max(temclust);

    center(b+1,1)=inde(ind);

end

for i=1:1:nopts
        disMin=inf;
        indexPointMin=0;
        for j=1:1:size(center,1)
            if Gdist(i,center(j,1))<disMin
                indexPointMin=j;
                disMin=Gdist(i,center(j,1));
            end
        end
        f_vec(1,i)=indexPointMin;
end

f_vec3=f_vec;

infile='centaur1.off';

fid=fopen(infile);
fgetl(fid);
nos=fscanf(fid, '%d %d %d',[3 1]);
nopts=nos(1);
notrg=nos(2);

coord=fscanf(fid, '%g %g %g',[3 nopts]);
coord=coord';

triang=fscanf(fid, '%d %d %d %d',[4 notrg]);
triang=triang';
triang=triang(:,2:4)+1;

fclose(fid);

    
outfile='centaur1kcenters.vtk';
ofid=fopen(outfile,'w');
fprintf(ofid, '# vtk DataFile Version 3.0\n');
fprintf(ofid, 'vtk output\n');
fprintf(ofid,'ASCII\n');
fprintf(ofid,'DATASET POLYDATA\n');
fprintf(ofid,'POINTS %d float\n',nopts);
fprintf(ofid,'%g %g %g\n',coord');
fprintf(ofid,'POLYGONS %d %d\n',notrg,4*notrg);
fprintf(ofid,'3 %d %d %d\n',triang'-1);
fprintf(ofid,'\n');
fprintf(ofid,'POINT_DATA %d\n',nopts);
fprintf(ofid,'SCALARS distance_from float 1\n');
fprintf(ofid,'LOOKUP_TABLE default\n');
fprintf(ofid,'%g\n',f_vec1);
fclose(ofid);

infile='centaur2.off';

fid=fopen(infile);
fgetl(fid);
nos=fscanf(fid, '%d %d %d',[3 1]);
nopts=nos(1);
notrg=nos(2);

coord=fscanf(fid, '%g %g %g',[3 nopts]);
coord=coord';

triang=fscanf(fid, '%d %d %d %d',[4 notrg]);
triang=triang';
triang=triang(:,2:4)+1;

fclose(fid);

    
outfile='centaur2kcenters.vtk';
ofid=fopen(outfile,'w');
fprintf(ofid, '# vtk DataFile Version 3.0\n');
fprintf(ofid, 'vtk output\n');
fprintf(ofid,'ASCII\n');
fprintf(ofid,'DATASET POLYDATA\n');
fprintf(ofid,'POINTS %d float\n',nopts);
fprintf(ofid,'%g %g %g\n',coord');
fprintf(ofid,'POLYGONS %d %d\n',notrg,4*notrg);
fprintf(ofid,'3 %d %d %d\n',triang'-1);
fprintf(ofid,'\n');
fprintf(ofid,'POINT_DATA %d\n',nopts);
fprintf(ofid,'SCALARS distance_from float 1\n');
fprintf(ofid,'LOOKUP_TABLE default\n');
fprintf(ofid,'%g\n',f_vec2);
fclose(ofid);


infile='centaur3.off';

fid=fopen(infile);
fgetl(fid);
nos=fscanf(fid, '%d %d %d',[3 1]);
nopts=nos(1);
notrg=nos(2);

coord=fscanf(fid, '%g %g %g',[3 nopts]);
coord=coord';

triang=fscanf(fid, '%d %d %d %d',[4 notrg]);
triang=triang';
triang=triang(:,2:4)+1;

fclose(fid);

    
outfile='centaur3kcenters.vtk';
ofid=fopen(outfile,'w');
fprintf(ofid, '# vtk DataFile Version 3.0\n');
fprintf(ofid, 'vtk output\n');
fprintf(ofid,'ASCII\n');
fprintf(ofid,'DATASET POLYDATA\n');
fprintf(ofid,'POINTS %d float\n',nopts);
fprintf(ofid,'%g %g %g\n',coord');
fprintf(ofid,'POLYGONS %d %d\n',notrg,4*notrg);
fprintf(ofid,'3 %d %d %d\n',triang'-1);
fprintf(ofid,'\n');
fprintf(ofid,'POINT_DATA %d\n',nopts);
fprintf(ofid,'SCALARS distance_from float 1\n');
fprintf(ofid,'LOOKUP_TABLE default\n');
fprintf(ofid,'%g\n',f_vec3);
fclose(ofid);