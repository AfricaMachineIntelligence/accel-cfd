function [indx] = calc_index_3d(i,j,k,m,n,l)

% Create Abbriviation for interior lengths
im = m-2;
in = n-2;
il = l-2;

% Create Abbrivation for blocks
b_clc = 2*4+im*5;
b_lsl = 2*5+im*6;
b_sis = 2*6+im*7;
plane1 = 2*b_clc+in*b_lsl;
plane2 = 2*b_lsl+in*b_sis;

% Run Index Check and Calculate
if k == 1
    if j == 1
        if i == 1
            indx = 1;  
        else
            indx = 5*(i-2)+4+2;
        end
    elseif j < n
        if i == 1
            indx = 2+(j-2)*b_lsl+b_clc;
        else
            indx = 6*(i-2)+8+(j-2)*b_lsl+b_clc;
        end
    else
        if i == 1
            indx = 2+in*b_lsl+b_clc;
        else
            indx = 5*(i-2)+7+in*b_lsl+b_clc;
        end
    end
    
elseif k < l
    if j == 1
        if i == 1
            indx = 2+(k-2)*plane2+plane1;  
        else
            indx = 6*(i-2)+8+(k-2)*plane2+plane1;
        end
    elseif j < n
        if i == 1
            indx = 3+(j-2)*b_sis+b_lsl+(k-2)*plane2+plane1;
        else
            indx = 7*(i-2)+10+(j-2)*b_sis+b_lsl+(k-2)*plane2+plane1;
        end
    else
        if i == 1
            indx = 3+in*b_sis+b_lsl+(k-2)*plane2+plane1;
        else
            indx = 6*(i-2)+9+in*b_sis+b_lsl+(k-2)*plane2+plane1;
        end
    end
    
else
    if j == 1
        if i == 1
            indx = 2+il*plane2+plane1;  
        else
            indx = 5*(i-2)+7+il*plane2+plane1;
        end
    elseif j < n
        if i == 1
            indx = 3+(j-2)*b_lsl+b_clc+il*plane2+plane1;
        else
            indx = 6*(i-2)+9+(j-2)*b_lsl+b_clc+il*plane2+plane1;
        end
    else
        if i == 1
            indx = 3+in*b_lsl+b_clc+il*plane2+plane1;
        else
            indx = 5*(i-2)+8+in*b_lsl+b_clc+il*plane2+plane1;
        end
    end
end
