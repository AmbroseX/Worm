function [Ydistance,Stopage] = calculate_mean_y_distance_N_stopage(individual_edge,bodyFraction_thrsh,min_y_thrshd,max_abs_y_distance)
coordinates = find_same_value(individual_edge); % index is cell, stores indexes have same each X.
if ~isempty(coordinates{1})    % Only extracts overlaped lineds.
    [Brims,Bnum,overlap_count] = if_contineous_line(coordinates,min_y_thrshd);  % Brims is cell, stores x,y coordinates of each Brim.
    if overlap_count >= bodyFraction_thrsh
        Bnumu = flip(unique(Bnum)); % flip the uniquing result, so get an array that sorts from largest to smallest. 
        Bnumu = Bnumu(Bnumu ~= 1);  % only operate overlaped points. 
        for j = 1:length(Bnumu)
            Brims_overlap = Brims(Bnum == Bnumu(j));
                if length(Brims_overlap) >= bodyFraction_thrsh
                    Brims_overlap = cell2mat(Brims_overlap);
                    distance = diff(Brims_overlap);
                    distance = distance(distance ~= 0);
                    Ydistance = mean(distance);
                    %Stopage = Brims_overlap(1,end-1)+3;
                    % find y limits for xStopage
                    jj = 0;
                    for ii = 1:length(Brims_overlap)
                        if ~mod(ii,2)
                            jj = jj+1;
                            yLimits(:,jj) = Brims_overlap(:,ii);
                        end
                    end
                    y_min = min(min(yLimits));
                    y_max = max(max(yLimits));
                    [y_individual_edge,x_individual_edge] = find(individual_edge);
                    xLimits = x_individual_edge(find((y_individual_edge <= y_max).*(y_min <= y_individual_edge)));
                    Stopage = max(xLimits);
                    if y_max - y_min > max_abs_y_distance
                        Stopage = 0;
                    end
                    break
                else
                    Ydistance = 0;
                    Stopage = 0;
                end
        end
    else
        Ydistance = 0;
        Stopage = 0;

    end
else
    Ydistance = 0;
    Stopage = 0;
end
end

function coordinates = find_same_value(individual_edge)
[y,x] = find(individual_edge);
xu = int16(unique(x));
coordinates = cell(1,1);
j = 0;
for i = 1:length(xu)
    index = x == xu(i);
    if length(find(index)) >1
    xcoordinates = x(index);
    ycoordinates = y(index);
    j = j+1;
    coordinates{j} = int16([xcoordinates,ycoordinates]);
    end
end
end
%% simply the line in single width and find seperate lines
function [Brims,Bnum,count] = if_contineous_line(coordinates,y_thrshd)
% y_threshd is the min gap to seperate two different brims.
for i = 1:length(coordinates)
    y = coordinates{i}(:,2);
    y_line = [0;y];
    % extract rows that are greater than y_threshold
    brimPoint = diff(y_line) > y_thrshd;
    Brims{i} = coordinates{i}(brimPoint,:);
    if y(1) <= y_thrshd;
        Brims{i} = [coordinates{i}(1,:);Brims{i}];
    end
end
count = 0; % number of overlaped points
for j = 1:length(Brims)
    Bnum(j) = size(Brims{j},1);
    if size(Brims{j},1) > 1
        count = count+1;
    end
end

end
    
        
            
