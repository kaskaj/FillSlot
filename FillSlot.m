function c = FillSlot(x,y,R,Sd,varname,var)
% INPUTS:
%           x,y - Coordinates of the slot polygon
%           R   - Conductor radius
%           Sd  - Density of conductor packing Sd is number between <0,1>
%                    Sd = 0 - Conductors are selected according to to weight w1
%                    Sd = 1 - Conductors are selected according to to weight w2
%                    w1 - Weight favoring the position closer to the bottom
%                    w2 - Weight favoring the position closer to center of the
%                         mass of the conductors
%           varname - Name of the parameter according to which the slot is filled
%                        varname = 'N'  - Fill slot with 'N' conductors
%                        varname = 'ff' - Fill slot with fill factor 'ff'
%                        (fill factor = area of conductors / area of slot)
%           var     - Parameter according to which the slot is filled
%                        if varname is 'N', var is the number of conductors 
%                        if varname is 'ff', var is the fill factor
% OUTPUTS:
%           c - Vector of x,y coordinates of conductor centers

%% Meta-Parameters

dh    = R/7;        % Mesh step (smaller step -> slower computation but more freedom)
t_SCI = 2*dh;       % Thickness of the initial strip of points
t_Rt  = 2*dh;       % Thickness of the strip of new points around new conductor
k     = 40;         % Level of curvature of the paraboloid defining the weight w2 

%% Check input parameters

if ~(strcmp(varname,'N') || strcmp(varname,'ff'))
    error('Unknown varname!');
end

if floor((max(y)-min(y))/dh) > 1000 || floor((max(x)-min(x))/dh) > 1000
    warning('The conductors are probably too small for the slot');
end  

%% Set regions

SC  = polyshape(x,y);                           % Whole space for conductors
SCR = polybuffer(SC,-R,'JointType','miter');    % Shrink SC by R
SCI = translate(SCR,0,-t_SCI);                  % Initial space for new center selection (default: bottom of the slot)

As = polyarea(SC.Vertices(:,1),SC.Vertices(:,2));      % Area of the slot
Ac = pi*R^2;                                           % Area of one conductor
 
%% Mesh for conductor centers

x_min = min(SCR.Vertices(:,1)); x_max = max(SCR.Vertices(:,1));
y_min_t = min(SCR.Vertices(:,2)); y_max_t = max(SCR.Vertices(:,2));

mesh = combvec(linspace(x_min,x_max,floor((x_max-x_min)/dh)), ...
               linspace(y_min_t,y_max_t,floor((y_max_t-y_min_t)/dh)))';
           
in = inpolygon(mesh(:,1),mesh(:,2),SCR.Vertices(:,1),SCR.Vertices(:,2));    % Points inside SCR
mesh(~in,:) = [];                                                           % Exclude points that do not belong to the SCR

%% Project regions on mesh

i_SCR =  inpolygon(mesh(:,1),mesh(:,2),SCR.Vertices(:,1),SCR.Vertices(:,2));     % Points of the mesh in SCR
i_SCI = ~inpolygon(mesh(:,1),mesh(:,2),SCI.Vertices(:,1),SCI.Vertices(:,2));     % Points of the mesh in SCI

Rn = ~i_SCR;                    % Rn - Vector of points from which a new center must not be selected
Rt =  i_SCR & i_SCI;            % Rt - Vector of points from which a new center can be selected

%% Loop settings

n    = 0;                   % Initial number of conductors
i_t  = find(Rt);            % List of free points
ff   = 0;                   % Initial fill factor
flag = 0;                   % Indication of premature end of loop

if strcmp(varname,'N'), i_c = zeros(var,1); else, i_c = zeros(floor(As/Ac),1); end      % Vector of conductor centers

%% Loop 

while ~(strcmp(varname,'N')  && n  == var || ...        % Stop condition for number of conductors
        strcmp(varname,'ff') && ff >= var)              % Stop condition for fill factor
    
    % If the stop conditions are not met, the loop stops after the entire slot is filled

   %% Generate new center
   
    n = n+1;                                        % Iterate
    
	if n == 1        
        i_c(n) = randsample(i_t,1);                 % Generate first center point
    else        
        r = rand(1);                                % Generate random number from < 0, 1 >         
        if r > Sd       
            if length(i_t) == 1, i_c(n) = i_t; else, i_c(n) = randsample(i_t,1,true,w1(i_t)); end     % Generate the center according to w1
        else
            if length(i_t) == 1, i_c(n) = i_t; else, i_c(n) = randsample(i_t,1,true,w2(i_t)); end     % Generate the center according to w2            
        end
    end
    
	%% Update regions
   
	Rn = Rn | (mesh(:,1) <  sqrt((2*R)^2 - (mesh(:,2) - mesh(i_c(n),2)).^2) + mesh(i_c(n),1) ...     
             & mesh(:,1) > -sqrt((2*R)^2 - (mesh(:,2) - mesh(i_c(n),2)).^2) + mesh(i_c(n),1));
	Rt = Rt | (mesh(:,1) <  sqrt((2*R + t_Rt)^2 - (mesh(:,2) - mesh(i_c(n),2)).^2) + mesh(i_c(n),1) ...     
             & mesh(:,1) > -sqrt((2*R + t_Rt)^2 - (mesh(:,2) - mesh(i_c(n),2)).^2) + mesh(i_c(n),1));
        
	Rt  = Rt & ~Rn;        % Update vectors Rn, Rt   
    ff  = n*(Ac/As);       % Update fill factor
	i_t = find(Rt);        % Update list of free points
    
    if isempty(i_t), flag = 1; break; end	% Break if there are no free points left
    
    %% Update weight w1 (prefer points closer to the bottom of the slot)
    
    y_max_t = max(mesh(Rt,2)); y_min_t = min(mesh(Rt,2));
    K1 = (1-0.1)/(y_max_t-y_min_t); K2 = 1-K1*y_max_t;      % Line w1 = K1*y + K2
    w1 = K1*mesh(:,2) + K2;                                 % Weight w1
    if isnan(w1(Rt)), w1(Rt) = 1; end                       % If y_max = y_min
    
    %% Update weight w2 (prefer points closer to the center of mass)
    
	mc  = sum(mesh(i_c(1:n),:),1)./n;       % Update center of mass 
	i_m = i_t(dsearchn(mesh(Rt,:),mc));     % Closest point of Rt to mc
 
    w2 = (-k/(R^2))*((mesh(:,1)-mesh(i_m,1)).^2 + (mesh(:,2)-mesh(i_m,2)).^2) + 1;       % Weight w2
    ii = w2 < 0; w2(ii) = 0;                                                             % Limit w2 to <0,1>

end

%% Coordinates of conductor centers

c  = mesh(i_c(1:n),:);       % Centers of conductors in (x,y) coordinates

%% Check flag

if flag, fprintf('The algorithm failed to reach the required number of conductors, try changing Sd.\n'); end
fprintf('Final number of conductors = %d, final fill factor = %d.\n', n, ff);

end