clear;

%% Dimensions of the trapezoid

b1 = 1;    % Bottom length
b2 = 2;    % Top length
h  = 3;    % Height

%% Radius of the conductor

R = 0.7e-1;

%% Trapezoid coordinates

x = [b1/2; b2/2; -b2/2; -b1/2];
y = [0; h; h; 0];

%% Fill the trapezoid

Sd      = 0.7;   % Density of conductor packing <0,1>
varname = 'N';   % varname = 'N'  - Fill slot with 'N' conductors
                 % varname = 'ff' - Fill slot with fill factor 'ff'
var     = 100;   % Number of conductors 'N' or fill factor 'ff'

c       = FillSlot(x,y,R,Sd,varname,var);

%% Plot result

PlotSlot(x,y,R,c);