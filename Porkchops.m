% error control
epsilon_newton = 10^-9;
epsilon_house1 = 10^-9;
epsilon_house2 = 10^-12;
epsilon_halley = 10^-12;
max_newton = 30;
max_house = 5;
max_halley = 10;


% system parameters
mu = 1.32712440018 * 10^20;
body1 = readmatrix('OrbitalElements.xlsx', 'Sheet', 'Sheet1', 'Range', 'B1:G21');
body1 = conversion2(body1);
body2 = readmatrix('OrbitalElements.xlsx', 'Sheet', 'Sheet1', 'Range', 'H1:M21');
body2 = conversion2(body2);
body3 = readmatrix('OrbitalElements.xlsx', 'Sheet', 'Sheet1', 'Range', 'N1:S21');
body3 = conversion2(body3);
body4 = readmatrix('OrbitalElements.xlsx', 'Sheet', 'Sheet1', 'Range', 'T1:Y21');
body4 = conversion2(body4);
body5 = readmatrix('OrbitalElements.xlsx', 'Sheet', 'Sheet1', 'Range', 'Z1:AE21');
body5 = conversion2(body5);
body6 = readmatrix('OrbitalElements.xlsx', 'Sheet', 'Sheet1', 'Range', 'AF1:AK21');
body6 = conversion2(body6);
body7 = readmatrix('OrbitalElements.xlsx', 'Sheet', 'Sheet1', 'Range', 'AL1:AQ21');
body7 = conversion2(body7);
body8 = readmatrix('OrbitalElements.xlsx', 'Sheet', 'Sheet1', 'Range', 'AR1:AW21');
body8 = conversion2(body8);
body9 = readmatrix('OrbitalElements.xlsx', 'Sheet', 'Sheet1', 'Range', 'AX1:BC21');
body9 = conversion2(body9);
body10 = readmatrix('OrbitalElements.xlsx', 'Sheet', 'Sheet1', 'Range', 'BD1:BI21');
body10 = conversion2(body10);
body11 = readmatrix('OrbitalElements.xlsx', 'Sheet', 'Sheet1', 'Range', 'BJ1:BO21');
body11 = conversion2(body11);
body12 = readmatrix('OrbitalElements.xlsx', 'Sheet', 'Sheet1', 'Range', 'BP1:BU21');
body12 = conversion2(body12);
bodymatrix = cat(3, body1, body2, body3, body4, body5, body6, body7, body8, body9, body10, body11, body12);


% day length increase per year due to tidal acceleration
DayLengthIncrease = 1.7 * 10^(-5);


% initial time for orbital elements matrix in Julian years from t = 0
tcelestial = 950;
maxjul = tcelestial + size(bodymatrix,1) - 1;
% spacing between each matrix row in Julian years
celestialinterval = 1;


% determines which trajectory to use, in case multiple valid trajectories exist
selectionmode = 1;
weight = 1;
% determines quantities shown on graph
displaymode = 6;
% determines rendering of graph
fillmode = 1;


% transfer between objects with the following IDs (object of origin first)
ID = [4 5];


% initial mission time for first possible launch (Gregorian time)
missiont0 = [964 0];
% spacing between each examined launch dates in days
writeinterval = 1;
% examined departure & arrival dates in days after t0 (must be a positive multiple of the write interval)
max_dep = 731;
max_arr = 1200;
axis1 = 0:(writeinterval):max_dep;
axis2 = 1:(writeinterval):max_arr;
% minimum time of flight examined in days (must be a positive multiple of the write interval)
tof_min = 30;
% maximum & minimum delta-v plotted
maxdv = 5 * 10^4;
mindv = 1;
% maximum & minimum C3 plotted
maxC3 = 500;
minC3 = 1;
lvls1 = 40;
lvls2 = 20;
% contour levels plotted
if or(((displaymode) == 5), ((displaymode) == 6))
    lvls1 = minC3:((maxC3)/(lvls1)):maxC3;
else
    lvls1 = mindv:((maxdv)/(lvls1)):maxdv;
end
if or(((displaymode) == 2), ((displaymode) == 5))
    lvls2 = mindv:((maxdv)/(lvls2)):maxdv;
end


maxcolumn = ((max_dep) / (writeinterval)) + 1;
maxrow = (max_arr) / (writeinterval);
out1 = zeros(maxrow, maxcolumn);
if or(((displaymode) == 2), ((displaymode) == 5))
    out2 = zeros(maxrow, maxcolumn);
end
columnindex = 1;
tod = 0;
while true
    t_dg = DTG(missiont0, tod);
    t_ds = GTE(t_dg);
    t_dj = ETJ(t_ds, maxjul);
    Kepler_dt = bodymatrix((t_dj(1) - tcelestial), :, (ID(1)));
    Kepler_d = transpose(Kepler_dt);
    S1 = KTC(Kepler_d, mu, t_dj(2), 0, epsilon_newton, max_newton);
    rowindex = (tod + tof_min) / (writeinterval);
    tof = tof_min;
    while true
        t_ag = DTG(t_dg, tof);
        t_as = GTE(t_ag);
        t_aj = ETJ(t_as, maxjul);
        Kepler_at = bodymatrix((t_aj(1) - tcelestial), :, (ID(2)));
        Kepler_a = transpose(Kepler_at);
        S2 = KTC(Kepler_a, mu, t_aj(2), 0, epsilon_newton, max_newton);
        deltat = t_as - t_ds;
        if or(((displaymode) == 2), ((displaymode) == 5))
            [output1, output2] = Lambert(S1, S2, deltat, mu, selectionmode, weight, displaymode, epsilon_house1, epsilon_house2, epsilon_halley, max_house, max_halley);
            if output1 > maxdv
                output1 = 0;
            end
            if displaymode == 5
                output1 = ((output1)^2) * 10^-6;
            end
            if output2 > maxdv
                output2 = 0;
            end
            out1(rowindex, columnindex) = output1;
            out2(rowindex, columnindex) = output2;
        else
            [output1, ~] = Lambert(S1, S2, deltat, mu, selectionmode, weight, displaymode, epsilon_house1, epsilon_house2, epsilon_halley, max_house, max_halley);
            if output1 > maxdv
                output1 = 0;
            end
            if displaymode == 6
                output1 = ((output1)^2) * 10^-6;
            end
            out1(rowindex, columnindex) = output1;
        end
        rowindex = rowindex + 1;
        if rowindex > maxrow
            break
        end
        tof = tof + writeinterval;
    end
    columnindex = columnindex + 1;
    if columnindex > maxcolumn
        break
    end
    tod = tod + writeinterval;
end


if fillmode == 1
    contourf(axis1, axis2, out1, lvls1)
    colormap(flipud(jet))
    shading interp
    colorbar
else
    contour(axis1, axis2, out1, lvls1)
end
if or(((displaymode) == 2), ((displaymode) == 5))
    hold on
    if fillmode == 2
        contourf(axis1, axis2, out2, lvls2)
        colormap(flipud(jet))
        shading interp
        colorbar
    else
        contour(axis1, axis2, out2, lvls2, "magenta", 'ShowText', 'on')
    end
end


function GreTime = DTG(GreTime, Plus)
GreTime = GreTime + [0 Plus];
if isinteger((GreTime(1)) / 100) == true
    MD = 365;
elseif isinteger((GreTime(1)) / 4) == true
    MD = 366;
else
    MD = 365;
end
if GreTime(2) > MD
    Cap = true;
    while Cap == true
        if isinteger((GreTime(1)) / 100) == true
            MD = 365;
        elseif isinteger((GreTime(1)) / 4) == true
            MD = 366;
        else
            MD = 365;
        end
        if GreTime(2) > MD
            GreTime = GreTime + [1 -MD];
        else
            Cap = false;
        end
    end
end
end


function SigTime = GTE(GreTime)
Yr900Time = (((86399.9813 + 86399.9997) / 2) * (292194)) + (((86399.9997 + 86400.002) / 2) * (36524));
Yrs = GreTime(1) - 900;
LYrs = floor(Yrs/4);
Yrs = Yrs - ((LYrs) * 4);
SigTime = (Yr900Time) + ((LYrs) * 4 * 365.25 * 86400.002) + ((Yrs) * 365 * 86400.002) + ((GreTime(2)) * 86400.002);
end


function JulTime = ETJ(SigTime, maxjul)
% finds Julian year and time from year start in seconds from time from t = 0
Year = floor((SigTime) / (365.25 * 86400));
Year = min(Year, maxjul);
Seconds = (SigTime) - ((Year) * (365.25 * 86400));
JulTime = [Year Seconds];
end


function newKepler = conversion1(Kepler)
    newKepler = Kepler;
    newKepler(1, :) = Kepler(1, :) * 149597870700;
    newKepler(3:6, :) = (pi / 180) * Kepler(3:6, :);
end


function newKepler = conversion2(Kepler)
    newKepler = Kepler;
    newKepler(:, 1) = Kepler(:, 1) * 149597870700;
    newKepler(:, 3:6) = (pi / 180) * Kepler(:, 3:6);
end


function S = KTC(Kepler, mu, t, t0, epsilon_newton, max_newton)
% convert Keplerian orbital elements stored as a vector into individual scalar quantities
% use negative values for a for hyperbolas
% use the value of the magnitude of the radius at periapsis for a for parabolas
a = Kepler(1);
e = Kepler(2);
i = Kepler(3);
long = Kepler(4);
arg = Kepler(5);
M0 = Kepler(6);

% determine current mean anomaly based on initial mean anomaly & elapsed time
if t == t0
    M = M0;
else
    deltat = t - t0;
    if e > 1
        M = M0 + deltat * sqrt(mu / (-(a)^3));
    elseif e ==1
        M = M0 + deltat * sqrt(mu / (2 * (a)^3));
    else
        M = M0 + deltat * sqrt(mu / ((a)^3));
    end
end

% find true anomaly
% for elliptical and hyperbolic orbits, call iterative Newton-Raphson method solver to iteratively solve for eccentric anomaly
if e > 1
    E_anomaly = Newton(e, M, epsilon_newton, max_newton);
    v_anomaly = 2 * atan(sqrt((e + 1) / (e - 1)) * tanh((E_anomaly) / 2));
elseif e == 1
    v_anomaly = 2 * atan(2 * sinh(asinh((3 / 2) * M) / 3));
else
    E_anomaly = Newton(e, M, epsilon_newton, max_newton);
    v_anomaly = 2 * atan(sqrt((1 + e) / (1 - e)) * tan((E_anomaly) / 2));
end

% determine semi-latus rectum
if e == 1
    l = 2 * a;
else
    l = a * (1 - (e)^2);
end

% determine the magnitude of the radius using semi-latus rectum, eccentricity, & true anomaly
rmag = l / (1 + e * cos(v_anomaly));

% use true anomaly to calculate unrotated radius vector in the orbital plane
xini = rmag * cos(v_anomaly);
yini = rmag * sin(v_anomaly);
rini = [xini; yini; 0];

% find perpendicular to the radius vector
rperp = cross([0; 0; 1], rini);
rperp = rperp / norm(rperp);

% find flight path angle
phi = (e * sin(v_anomaly)) / (1 + e * cos(v_anomaly));

% find magnitude of the velocity vector using the vis-visa equation
if e == 1
    vmag = sqrt((2 * mu) / rmag);
else
    vmag = sqrt(mu * ((2 / rmag) - (1 / a)));
end

% find direction of velocity vector in the orbital plane
vini = [cos(-phi) -sin(-phi) 0; sin(-phi) cos(-phi) 0; 0 0 1] * rperp;

% find unrotated velocity vector in the orbital plane
vini = vmag * vini;

% apply rotation in the orbital plane from the argument of periapsis
rrot = [cos(arg) -sin(arg) 0; sin(arg) cos(arg) 0; 0 0 1] * rini;
vrot = [cos(arg) -sin(arg) 0; sin(arg) cos(arg) 0; 0 0 1] * vini;

% find orbital plane orientation in the outside reference frame
worldmat = [cos(long) (-sin(long)*cos(i)) (sin(long)*sin(i)); sin(long) (cos(long)*cos(i)) -(cos(long)*sin(i)); 0 sin(i) cos(i)];

% find radius and velocity vectors in the outside reference frame
r = worldmat * rrot;
v = worldmat * vrot;

% find state vector
S = [r; v];
end


function E_anomaly = Newton(e, M, epsilon_newton, max_newton)
% set error threshold
epsilon = epsilon_newton;

% initialize first guess
E_temp = M;

% use Newton-Raphson method to find the roots of Kepler's equation
if e < 1
    for i = 1:max_newton
        E_new = E_temp - ((E_temp - (e * sin(E_temp)) - M) / (1 - (e * cos(E_temp))));
        if (E_new - E_temp) <= epsilon
            E_temp = E_new;
            break
        end
        E_temp = E_new;
    end
else
    for i = 1:max_newton
        E_new = E_temp - (((e * sinh(E_temp)) - E_temp - M) / ((e * cosh(E_temp)) - 1));
        if (E_new - E_temp) <= epsilon
            E_temp = E_new;
            break
        end
        E_temp = E_new;
    end
end
E_anomaly = E_temp;
end


function Kepler = CTK(S, mu)
r = S(1:3);
v = S(4:6);

% find magnitudes of the radius & velocity vectors
rmag = norm(r);
vmag = norm(v);

% find angular momentum vector & magnitude
h = cross(r, v);
hmag = norm(h);

% find semi-latus rectum
l = ((hmag)^2) / mu;

% find semi-major axis& eccentricity
if vmag == sqrt((2 * mu) / rmag)
    a = l / 2;
    e = 1;
else
    a = -(1 / ((((vmag)^2) / mu) - (2 / rmag)));
    e = sqrt(1 - (l / a));
end

% find longitude of the rising node
hx = h(1);
hy = h(2);
long = atan2((hx), -(hy));

% find inclination
p = -(cross(([0; 0; 1]), ([cos(long); sin(long); 0])));
zproj = dot(([0; 0; 1]), (h));
pproj = dot((p), (h));
i = atan2(pproj, zproj);

% find true anomaly & argument of periapsis
if and((i == 0), (e == 0))
    x = r(1);
    y = r(2);
    v_anomaly = atan2(y, x);
    arg = 0;
elseif e == 0
    x = r(1);
    y = r(2);
    z = r(3);
    v_anomaly = atan2((z / (sin(i))), ((x * (cos(long))) + (y * (sin(long)))));
    arg = 0;
else
    x = r(1);
    y = r(2);
    z = r(3);
    v_anomaly = atan2((sqrt(l / (mu)) * dot(v, r)), (l - rmag));
    arglat = atan2((z / (sin(i))), ((x * (cos(long))) + (y * (sin(long)))));
    arg = arglat - v_anomaly;
end

% find mean anomaly
if e > 1
    E_anomaly = 2 * atanh(sqrt((e - 1) / (e + 1)) * tan((v_anomaly) / 2));
    M = (e * sinh(E_anomaly)) - (E_anomaly);
elseif e == 1
    M = (2 / 3) * sinh(3 * asinh((tan((v_anomaly) / 2)) / 2));
else
    E_anomaly = 2 * atan(sqrt((1 - e) / (e + 1)) * tan((v_anomaly) / 2));
    M = (E_anomaly) - (e * sin(E_anomaly));
end

% assemble vector from orbital elements
Kepler = [a e i long arg M];
end


function [output1, output2] = Lambert(S1, S2, t, mu, selectionmode, weight, displaymode, epsilon_house1, epsilon_house2, epsilon_halley, max_house, max_halley)
arguments
    S1
    S2
    t {mustBeGreaterThan(t, 0)}
    mu {mustBeGreaterThan(mu, 0)}
    selectionmode {mustBeInteger, mustBeInRange(selectionmode, 1, 4)}
    weight
    displaymode {mustBeInteger, mustBeInRange(displaymode, 1, 6)}
    epsilon_house1
    epsilon_house2
    epsilon_halley
    max_house {mustBeInteger}
    max_halley {mustBeInteger}
end
% extract initial and final radius and velocity
r1 = S1(1:3);
v1 = S1(4:6);
r2 = S2(1:3);
v2 = S2(4:6);
c = r2 - r1;
cnorm = norm(c);
r1norm = norm(r1);
r2norm = norm(r2);
s = (1 / 2) * (r1norm + r2norm + cnorm);

% radial unit vectors
a1 = r1 / r1norm;
a2 = r2 / r2norm;

% orbital momentum and central angle
h = cross(a1, a2);
hnorm = norm(h);
h = h / hnorm;
lambdasquared = 1 - (cnorm / s);
lambda = sqrt(lambdasquared);

% tangential unit vectors
if ((r1(1) * r2(2)) - (r1(2) * r2(1))) < 0
    lambda = -lambda;
    b1 = cross(a1, h);
    b2 = cross(a2, h);
else
    b1 = cross(h, a1);
    b2 = cross(h, a2);
end

% define new time
T = sqrt((2 * mu) / ((s)^3)) * t;

% call householder iterative solver
SolutionsList = findxy(lambda, T, epsilon_house1, epsilon_house2, epsilon_halley, max_house, max_halley);

% find departure and arrival velocities
gamma = sqrt((mu * s) / 2);
rho = (r1norm - r2norm) / cnorm;
sigma = sqrt(1 - (rho)^2);
matrixsize = size(SolutionsList, 1);
VelocityDeparture = zeros(3, matrixsize);
VelocityArrival = zeros(3, matrixsize);
v1 = repmat(v1, 1, matrixsize);
v2 = repmat(v2, 1, matrixsize);
for i = 1:matrixsize
    x = SolutionsList(i, 1);
    y = SolutionsList(i, 2);
    vr1 = gamma * (((lambda * y) - x) - (rho * ((lambda * y) + x))) / r1norm;
    vr2 = -gamma * (((lambda * y) - x) + (rho * ((lambda * y) + x))) / r2norm;
    vt1 = gamma * sigma * (y + (lambda * x)) / r1norm;
    vt2 = gamma * sigma * (y + (lambda * x)) / r2norm;
    vdep = vr1 * a1 + vt1 * b1;
    varr = vr2 * a2 + vt2 * b2;
    VelocityDeparture(:, i) = vdep;
    VelocityArrival(:, i) = varr;
end

% subtract initial and final velocities
RelVelDep = VelocityDeparture - v1;
RelVelArr = VelocityArrival - v2;

% if multiple solutions exist, generate costmatrix
costmatrix = zeros(1, matrixsize);
if selectionmode == 1
    for i = 1:matrixsize
        costmatrix(1, i) = norm(RelVelDep(:, i)) + norm(RelVelArr(:, i));
    end
elseif selectionmode == 2
    for i = 1:matrixsize
        costmatrix(1, i) = weight * norm(RelVelDep(:, i)) + (1 - weight) * norm(RelVelArr(:, i));
    end
elseif selectionmode == 3
    for i = 1:matrixsize
        costmatrix(1, i) = (norm(RelVelDep(:, 1)))^2 + (norm(RelVelArr(:, 1)))^2;
    end
else
    for i = 1:matrixsize
        costmatrix(1, i) = weight * (norm(RelVelDep(:, 1)))^2 + (1 - weight) * (norm(RelVelArr(:, 1)))^2;
    end
end

% index optimal trajectory in costmatrix based on user-defined criteria
[~, Index] = min(costmatrix,[],"linear");
Index = Index(1);

% output user-defined cost of trajectory
if displaymode == 1
    output1 = norm(RelVelDep(:, Index)) + norm(RelVelArr(:, Index));
    output2 = 0;
elseif or((displaymode == 2), (displaymode == 5))
    output1 = norm(RelVelDep(:, Index));
    output2 = norm(RelVelArr(:, Index));
elseif or((displaymode == 3), (displaymode == 6))
    output1 = norm(RelVelDep(:, Index));
    output2 = 0;
elseif displaymode == 4
    output1 = norm(RelVelArr(:, Index));
    output2 = 0;
end
end


function SolutionsList = findxy(lambda, T, epsilon_house1, epsilon_house2, epsilon_halley, max_house, max_halley)
arguments
    lambda {mustBeLessThan(lambda, 1), mustBeGreaterThan(lambda, -1)}
    T {mustBeGreaterThan(T, 0)}
    epsilon_house1
    epsilon_house2
    epsilon_halley
    max_house
    max_halley
end
% maximum number of orbits before intercept
M_max = floor(T / pi);

% time of flight for x = 0 for zero-orbit intercept
T00 = acos(lambda) + lambda * sqrt(1 - (lambda)^2);
if T < T00 + M_max * pi && M_max > 0
    % call iterative Halley solver for minimum time-of-flight for maximum number of orbits intercept
    Tmin = Halley(0, M_max, lambda, T00, epsilon_halley, max_halley);
    if Tmin > T
        M_max = M_max - 1;
    end
end

% parabolic time of flight
T1 = (2/3) * (1 - (lambda)^3);
if T == T1
    x = 1;
    y = 1;
elseif T == T00
    x = 0;
    y = sqrt(1 - (lambda)^2);
else
    % initialize iterative solver for zero-orbit intercept
    if T >= T00
        x0 = (T00 / T)^(2/3) - 1;
    elseif T < T1
        x0 = (5/2) * ((T1 * (T1 - T)) / (T * (1 - (lambda)^5))) + 1;
    else
        x0 = (T00 / T)^(log2(T1 / T00)) - 1;
    end
    
    % call iterative Householder solver
    x = Householder(x0, 0, lambda, T, epsilon_house1, epsilon_house2, max_house);
    y = sqrt(1 - ((lambda)^2) * (1 - (x)^2));
end
if M_max > 0
    SolutionSize = 1 + 2 * (M_max);
    SolutionsList = zeros(SolutionSize,2);
end
SolutionsList(1, :) = [x y];
if M_max > 0
    for i = 1:M_max
        % initialize iterative solver for multi-orbit left-branch intercept
        M = i;
        x0l = ((((M * pi) + pi) / (8 * T))^(2 / 3) - 1) / ((((M * pi) + pi) / (8 * T))^(2 / 3) + 1);

        % call iterative Householder solver for left-branch intercept
        xl = Householder(x0l, M, lambda, T, epsilon_house1, epsilon_house2, max_house);
        yl = sqrt(1 - ((lambda)^2) * (1 - (xl)^2));
        SolutionsList((M * 2), :) = [xl yl];

        % initialize iterative solver for multi-orbit right-branch intercept
        x0r = (((8 * T) / (M * pi))^(2/3) - 1) / (((8 * T) / (M * pi))^(2/3) + 1);

        % call iterative Householder solver for right-branch intercept
        xr = Householder(x0r, M, lambda, T, epsilon_house1, epsilon_house2, max_house);
        yr = sqrt(1 - ((lambda)^2) * (1 - (xr)^2));
        SolutionsList(((M * 2) + 1), :) = [xr yr];
    end
end
end


function x = Householder(x0, M, lambda, Tactual, epsilon_house1, epsilon_house2, max_house)
% set error threshold
if M > 0
    epsilon = epsilon_house2;
else
    epsilon = epsilon_house1;
end
x = x0;

% use Householder's (third order) method to solve for trajectory
for i = 1:max_house
    y = sqrt(1 - ((lambda)^2) * (1 - (x)^2));
    if x < 1
        psi = acos((x * y) + lambda * (1 - (x)^2));
    else
        psi = acosh((x * y) - lambda * ((x)^2 - 1));
    end
    T = (1 / (1 - (x)^2)) * (((psi + M * pi) / (sqrt(abs(1 - (x)^2)))) - x + (lambda) * y);
    Tprime = ((3 * T * x) - 2 + (2 * (lambda)^3 * (x / y))) / (1 - (x)^2);
    Tdoubleprime = ((3 * T) + (5 * x * (Tprime)) + (2 * (1 - (lambda)^2) * ((lambda)^3 / (y)^3))) / (1 - (x)^2);
    Ttripleprime = ((7 * x * (Tdoubleprime)) + (8 * Tprime) - (6 * (1 - (lambda)^2) * (lambda)^5 * (x / ((y)^5)))) / (1 - (x)^2);
    Objective = (T - Tactual);
    xnew = x - Objective * (((Tprime)^2 - (((Objective) * (Tdoubleprime)) / 2)) / ((Tprime * (((Tprime)^2) - ((Objective) * (Tdoubleprime)))) + (((Ttripleprime) * (Objective)^2) / 6)));
    if (xnew - x) <= epsilon
        x = xnew;
        break
    end
    x = xnew;
end
end


function Tmin = Halley(x0, M, lambda, T00, epsilon_halley, max_halley)
% set error threshold
epsilon = epsilon_halley;
x = x0;

% use Halley's method to solve for minimum time-of-flight
for i = 1:max_halley
    y = sqrt(1 - ((lambda)^2) * (1 - (x)^2));
    psi = acos((x * y) + lambda * (1 - (x)^2));
    if i == 1
        T = T00 + (M * pi);
    else
        T = (1 / (1 - (x)^2)) * (((psi + M * pi) / (sqrt(abs(1 - (x)^2)))) - x + (lambda) * y);
    end
    Tprime = ((3 * T * x) - 2 + (2 * (lambda)^3 * (x / y))) / (1 - (x)^2);
    Tdoubleprime = ((3 * T) + (5 * x * (Tprime)) + (2 * (1 - (lambda)^2) * ((lambda)^3 / (y)^3))) / (1 - (x)^2);
    Ttripleprime = ((7 * x * (Tdoubleprime)) + (8 * Tprime) - (6 * (1 - (lambda)^2) * (lambda)^5 * (x / ((y)^5)))) / (1 - (x)^2);
    xnew = x - ((2 * (Tprime) * (Tdoubleprime)) / ((2 * (Tdoubleprime)^2) - ((Tprime) * (Ttripleprime))));
    if (xnew - x) <= epsilon
        x = xnew;
        break
    end
    x = xnew;
end
y = sqrt(1 - ((lambda)^2) * (1 - (x)^2));
psi = acos((x * y) + lambda * (1 - (x)^2));
Tmin = (1 / (1 - (x)^2)) * (((psi + M * pi) / (sqrt(abs(1 - (x)^2)))) - x + (lambda) * y);
end