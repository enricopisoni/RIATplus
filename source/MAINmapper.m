function []=MAINmapper(aqi,outdir,pathPRB,pathPCI,livelli)

% -aqi:     the map you want to create (column1=xutm column2=yutm col3=values)
% -outdir:  output directory and filename (i.e. outdir='./output/TestMaps/meanPM10';)
% -pathPRB: path for regional border (i.e. pathPRB='./inputOPERAaspa_LUG2013/regionsBoundaries.bln';)
% -pathPCI: path for cities (i.e. pathPCI='./inputOPERAaspa_LUG2013/citta.dat';)
% -livelli: colorbar levels (es 0:5:45 per PM10 oppure 10000:5000:90000 per AOT40)


figure;
%UNIBS(ET)20131002 - to transform coordinates are in m if they
%are in Km
if aqi(1,1)>1000
    unit=1;
else
    unit=1000;
end

POSTcomputeMapsRead(livelli,aqi,unit)
hold on

%draw regional boundaries
POSTcomputeMapsBln(pathPRB);

%draw cities
POSTcomputeMapsPost(pathPCI);
set(gcf,'Position',[1 31 1258 694]);
set(gcf,'PaperPositionMode','auto');

%save map of results
file=strcat(outdir,'.png');
print(gcf,'-dpng',file);
close

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FUNCTION ANCILLARY FOR MAPS CREATION
function POSTcomputeMapsRead(livelli,dati,units,~)
%contours of optimal aqi and emissions

%create coordinates
dx=dati(1,1);
diversi=find(dati(:,1)-dx>0.001);
%         diversi=find(dati(:,1)~=dx); %OLD RIAT
dx=round(abs(dx-dati(diversi(1),1)));
%         dx=abs(dx-dati(diversi(1),1));

dy=dati(1,2);
diversi=find(dati(:,2)-dy>0.001);
%         diversi=find(dati(:,2)~=dy);
dy=round(abs(dy-dati(diversi(1),2)));
%         dy=abs(dy-dati(diversi(1),2));

xorg=min(dati(:,1));
yorg=min(dati(:,2));

[nr nc]=size(dati);

%create map
for i=1:nr
    ii=round((dati(i,1)-xorg)/dx+1);
    jj=round((dati(i,2)-yorg)/dy+1);
    %             ii=(dati(i,1)-xorg)/dx+1; %OLD RIAT
    %             jj=(dati(i,2)-yorg)/dy+1;
    xcoord(ii)=(dati(i,1)+dx/2)*units;
    ycoord(jj)=(dati(i,2)+dy/2)*units;
    mappa(jj,ii)=dati(i,3);
end

%create contour
if isempty(livelli)==1
    contourf(xcoord,ycoord,mappa);
else
    contourf(xcoord,ycoord,mappa,livelli);
end
set(gca,'fontsize',18)
shading flat

if isempty(livelli)==1
    h=contourcmap('jet','colorbar','on');
else
    h=contourcmap(livelli,'jet','colorbar','on');
end

set(h,'fontsize',18)

y_tick=get(h,'ytick');
y_tick=y_tick(1:2:length(y_tick));
y_tick_label=get(h,'yticklabel');
y_tick_label=y_tick_label(1:2:length(y_tick_label),:);
set(h,'ytick',y_tick,'yticklabel',y_tick_label);

minimo=min(dati(:,1:2));
massimo=max(dati(:,1:2));

set(gca,'XLim',[minimo(1)*units massimo(1)*units]);
set(gca,'YLim',[minimo(2)*units massimo(2)*units]);

asp_ratio=get(gca,'DataAspectRatio');
asp_ratio(2)=asp_ratio(1);
set(gca,'DataAspectRatio',asp_ratio);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FUNCTION ANCILLARY FOR MAPS CREATION
function POSTcomputeMapsBln(infile)
%depict domain and regional boundaries

% read and plot bnl file
hold on
dati=load (infile);
blocchi_dati=find(dati(:,2)==0);
inizio=0;

for i=1:length(blocchi_dati)
    dati_plot=dati(2+inizio:2+inizio+dati(blocchi_dati(i),1)-1,:);
    inizio=inizio+dati(blocchi_dati(i),1)+1;
    plot(dati_plot(:,1),dati_plot(:,2),'Color',[0.5 0.5 0.5])
end

dati(blocchi_dati,:)=[];

minimo=min(dati);
massimo=max(dati);

asp_ratio=get(gca,'DataAspectRatio');
asp_ratio(2)=asp_ratio(1);
set(gca,'DataAspectRatio',asp_ratio);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FUNCTION ANCILLARY FOR MAPS CREATION
function POSTcomputeMapsPost(infile)
%depict cities

%open file
fid=fopen(infile);

while (feof(fid)==0)
    riga=fgetl(fid);
    if (feof(fid)==0)
        coord_xy(1,:)=sscanf(riga,'%f');
        blank=regexpi(riga,'\s');
        nome_label{1}=strtrim(riga(blank(length(blank)):length(riga)));
        
        h=text(coord_xy(:,1),coord_xy(:,2),strcat('+ ',nome_label{1}));
        set(h,'Fontsize',12);
    end
end
end

