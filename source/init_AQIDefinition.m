% inizializza le informazioni relative alla rete neurale
% patANN specializzato per contenere informazioni!
function [pathANN, pathDd, pathBc]=init_AQIDefinition(pathDEF, AQINum, aggType)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       aqi_definition.txt

% preparing the data structure to store the AQIs configuration
% 7 columns are the 7 AQIs
% 3 structures are: 1)annual, 2)winter, 3)summer

pathANN = struct ('ANNs',{});
pathDd = struct ('Dd',{});
pathBc = struct ('Bc',{});
for h=1:3;
    path = cellstr('-999');
    for p=1:AQINum-1,
        path = [path, cellstr('-999')];
    end;
    path = char(path);
    
    pathANN = [pathANN, struct('ANNs',path)];
    pathDd = [pathDd, struct('Dd',path)];
    pathBc = [pathBc, struct('Bc',path)];
end;

fid = fopen(pathDEF);
C1 = textscan(fid,'%s','delimiter','\n');

row = 1;
entryNum = str2num(cell2mat(C1{1}(row))); row = row + 1;

row = row + 1;
for k=1:entryNum,
    id = regexp(cell2mat(C1{1}(row)),'^[0-9]+','match');   row = row + 1;
    months = regexp(cell2mat(C1{1}(row)),'^[0-9]+','match');   row = row + 1;
    %     ann = regexp(cell2mat(C1{1}(row)),'\.*(/\w+)+','match');  row = row + 1;
    %     Dd = regexp(cell2mat(C1{1}(row)),'\.*(/\w+)+','match');  row = row + 1;
    %20130731 - generalize the use of paths - both absolute and relative
    ann=sscanf(cell2mat(C1{1}(row)),'%s %'); row = row + 1;
    Bc='-999';
    if (isequal(strtrim(aggType),'FIRSTGUESS')==1)
        Bc=sscanf(cell2mat(C1{1}(row)),'%s %');
        row = row + 1;
    end
    Dd=sscanf(cell2mat(C1{1}(row)),'%s %'); row = row + 1;
    row = row + 1;
    
    h = 0;
    if (isequal(strtrim(char(months)),'1')==1)
        h = 1;
    end
    if (isequal(strtrim(char(months)),'2')==1)
        h = 2;
    end
    if (isequal(strtrim(char(months)),'3')==1)
        h = 3;
    end
    
    jj = cell2mat(id)-char('0') + 1;
    if (jj == 1)
        pathANN(h) = struct('ANNs', char(char(ann), pathANN(h).ANNs(2:AQINum,:)));
        pathDd(h) = struct('Dd', char(char(Dd), pathDd(h).Dd(2:AQINum,:)));
        if (isequal(strtrim(aggType),'FIRSTGUESS')==1)
            pathBc(h)=struct('Bc', char(char(Bc), pathBc(h).Bc(2:AQINum,:)));
        end
    end
    if(jj == AQINum)
        pathANN(h) = struct('ANNs', char(pathANN(h).ANNs(1:AQINum-1,:), char(ann)));
        pathDd(h) = struct('Dd', char(pathDd(h).Dd(1:AQINum-1,:), char(Dd)));
        if (isequal(strtrim(aggType),'FIRSTGUESS')==1)
            pathBc(h) = struct('Bc', char(pathBc(h).Bc(1:AQINum-1,:), char(Bc)));
        end
    end
    if((jj > 1) && (jj < AQINum))
        pathANN(h) = ...
            struct('ANNs', char(pathANN(h).ANNs(1:jj-1,:), ...
            char(ann), ...
            pathANN(h).ANNs(jj+1:AQINum,:)));
        pathDd(h) = ...
            struct('Dd', char(pathDd(h).Dd(1:jj-1,:), ...
            char(Dd), ...
            pathDd(h).Dd(jj+1:AQINum,:)));
        if (isequal(strtrim(aggType),'FIRSTGUESS')==1)
            pathBc(h)= ...
                struct('Bc', char(pathBc(h).Bc(1:jj-1,:), ...
                char(Bc), ...
                pathBc(h).Bc(jj+1:AQINum,:)));
        end
    end
end;

fclose(fid);
end
