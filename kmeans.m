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

centaur1_hks=horzcat(centaur1_hks,coord);



for g=1:1:50
    center=randi(nopts,10,1);
    oldcenter=zeros(10,1);
while(isequal(center,oldcenter)==0)
oldcenter=center;

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
    
for k = 1:1:10
	ind1=find(f_vec==k);
    temp_hks=[];
    for i=1:1:size(ind1,2)
        temp_hks(i,:)=centaur1_hks(ind1(1,i),:);
    end
    temclust=mean(temp_hks);
    centmean=Euclidean_distance(temclust,temp_hks);
    [A i]=min(centmean);
    if(isempty(ind1)==0)
    center(k,1)=ind1(1,i);
    for i=1:1:size(ind1,2)
        disto(1,ind1(1,i))=Euclidean_distance(centaur1_hks(ind1(1,i),:),centaur1_hks(center(k,1),:));
    end
    end
end
end
distort(1,g)=sum(disto)/3400;
fvecl(g,:)=f_vec(1,:);
end

[a,i]=min(distort);
f_vec(1,:)=fvecl(i,:);
f_vec_1=f_vec;

fprintf('Program paused. Press enter to apply k center to the centaur 2\n');
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

centaur2_hks=horzcat(centaur2_hks,coord);



for g=1:1:50
    center=randi(nopts,10,1);
    oldcenter=zeros(10,1);
while(isequal(center,oldcenter)==0)
oldcenter=center;

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
    
for k = 1:1:10
	ind1=find(f_vec==k);
    temp_hks=[];
    for i=1:1:size(ind1,2)
        temp_hks(i,:)=centaur2_hks(ind1(1,i),:);
    end
    temclust=mean(temp_hks);
    centmean=Euclidean_distance(temclust,temp_hks);
    [A i]=min(centmean);
    if(isempty(ind1)==0)
    center(k,1)=ind1(1,i);
    for i=1:1:size(ind1,2)
        disto(1,ind1(1,i))=Euclidean_distance(centaur2_hks(ind1(1,i),:),centaur2_hks(center(k,1),:));
    end
    end
end
end
distort(1,g)=sum(disto)/3400;
fvecl(g,:)=f_vec(1,:);
end

[a,i]=min(distort);
f_vec(1,:)=fvecl(i,:);
f_vec_2=f_vec;

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

centaur3_hks=horzcat(centaur3_hks,coord);



for g=1:1:50
    center=randi(nopts,10,1);
    oldcenter=zeros(10,1);
while(isequal(center,oldcenter)==0)
oldcenter=center;

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
    
for k = 1:1:10
	ind1=find(f_vec==k);
    temp_hks=[];
    for i=1:1:size(ind1,2)
        temp_hks(i,:)=centaur3_hks(ind1(1,i),:);
    end
    temclust=mean(temp_hks);
    centmean=Euclidean_distance(temclust,temp_hks);
    [A i]=min(centmean);
    if(isempty(ind1)==0)
    center(k,1)=ind1(1,i);
    for i=1:1:size(ind1,2)
        disto(1,ind1(1,i))=Euclidean_distance(centaur3_hks(ind1(1,i),:),centaur3_hks(center(k,1),:));
    end
    end
end
end
distort(1,g)=sum(disto)/3400;
fvecl(g,:)=f_vec(1,:);
end

[a,i]=min(distort);
f_vec(1,:)=fvecl(i,:);
f_vec_3=f_vec;

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

    
outfile='centaur1kmeans.vtk';
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
fprintf(ofid,'%g\n',f_vec_1);
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

    
outfile='centaur2kmeans.vtk';
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
fprintf(ofid,'%g\n',f_vec_2);
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

    
outfile='centaur3kmeans.vtk';
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
fprintf(ofid,'%g\n',f_vec_3);
fclose(ofid);
