function [quad]=interface_get_quadrant(aggregationInfo)

 if isequal(aggregationInfo.type, 'FIRSTGUESS')
     quad =1;
 else
     % set to 4 with Quadrant
     quad=4;
 end

end
