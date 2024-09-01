% error control
epsilon_newton = 10^-9;
epsilon_house1 = 10^-9;
epsilon_house2 = 10^-12;
epsilon_halley = 10^-12;
epsilon_geoalt = 10^-12;
max_newton = 30;
max_house = 5;
max_halley = 10;
max_geoalt = 30;


% system parameters
mu = 1.075938*10^23*6.67408*10^-11;
a_body = 2.4103*10^6;                                         % semi-major axis (metres)
f_inverse = 298.257;                                        % inverse flattening
f_body = 0;
b_body = a_body * (1 - f_body);
e_first = sqrt(1 - ((b_body)^2 / (a_body)^2));

% sidereal rotation rate from Dehant & Matthews
omega_0 = (2 * pi) / (16.6890184 * 86400);

% secular precession and obliquity terms taken from Laskar's NGT
% secular precession terms
psi_0   = 0;
psi_1   = ((502909.66 * (1 / 3600)) * (pi / 180));
psi_2   = ((11119.71 * (1 / 3600)) * (pi / 180));
psi_3   = ((77.32 * (1 / 3600)) * (pi / 180));
psi_4   = (((-2353.16) * (1 / 3600)) * (pi / 180));
psi_5   = (((-180.55) * (1 / 3600)) * (pi / 180));
psi_6   = ((174.51 * (1 / 3600)) * (pi / 180));
psi_7   = ((130.95 * (1 / 3600)) * (pi / 180));
psi_8   = ((24.24 * (1 / 3600)) * (pi / 180));
psi_9   = (((-47.59) * (1 / 3600)) * (pi / 180));
psi_10  = ((-8.66 * (1 / 3600)) * (pi / 180));
psi_vec = [psi_0 psi_1 psi_2 psi_3 psi_4 psi_5 psi_6 psi_7 psi_8 psi_9 psi_10];
% secular obliquity terms
% (23 + (26 * (1 / 60)) + (21.448 * (1 / 3600))) * (pi / 180);
theta_0   = 0;
theta_1   = (((-4680.93) * (1 / 3600)) * (pi / 180));
theta_2   = (((-1.55) * (1 / 3600)) * (pi / 180));
theta_3   = ((1999.25 * (1 / 3600)) * (pi / 180));
theta_4   = (((-51.38) * (1 / 3600)) * (pi / 180));
theta_5   = (((-249.67) * (1 / 3600)) * (pi / 180));
theta_6   = (((-39.05) * (1 / 3600)) * (pi / 180));
theta_7   = ((7.12 * (1 / 3600)) * (pi / 180));
theta_8   = ((27.87 * (1 / 3600)) * (pi / 180));
theta_9   = ((5.79 * (1 / 3600)) * (pi / 180));
theta_10  = ((2.45 * (1 / 3600)) * (pi / 180));
theta_vec = [theta_0 theta_1 theta_2 theta_3 theta_4 theta_5 theta_6 theta_7 theta_8 theta_9 theta_10];
sec_matrix = [psi_vec; theta_vec];


% spectral precession terms from Dehant & Matthews
psi_ip1 = 0;
psi_ip2 = 0;
psi_ip3 = 0;
psi_ip4 = 0;
psi_ip5 = 0;
psi_ipvec = [psi_ip1; psi_ip2; psi_ip3; psi_ip4; psi_ip5];
psi_oop1 = 0;
psi_oop2 = 0;
psi_oop3 = 0;
psi_oop4 = 0;
psi_oop5 = 0;
psi_oopvec = [psi_oop1; psi_oop2; psi_oop3; psi_oop4; psi_oop5];
% spectral obliquity terms from Dehant & Matthews
theta_ip1 = 0;
theta_ip2 = 0;
theta_ip3 = 0;
theta_ip4 = 0;
theta_ip5 = 0;
theta_ipvec = [theta_ip1; theta_ip2; theta_ip3; theta_ip4; theta_ip5];
theta_oop1 = 0;
theta_oop2 = 0;
theta_oop3 = 0;
theta_oop4 = 0;
theta_oop5 = 0;
theta_oopvec = [theta_oop1; theta_oop2; theta_oop3; theta_oop4; theta_oop5];
% frequency of spectral nutation from Dehant & Matthews
omega_1 = 0;
omega_2 = 0;
omega_3 = 0;
omega_4 = 0;
omega_5 = 0;
omega_vec = [omega_1; omega_2; omega_3; omega_4; omega_5];
% phase component of spectral nutation from Dehant & Matthews
chi_1 = 0;
chi_2 = 0;
chi_3 = 0;
chi_4 = 0;
chi_5 = 0;
chi_vec = [chi_1; chi_2; chi_3; chi_4; chi_5];
spec_matrix = [psi_ipvec psi_oopvec theta_ipvec theta_oopvec omega_vec chi_vec];


% determines which trajectory to use, in case multiple valid trajectories exist
selectionmode = 1;
weight = 1;
% determines quantities shown on graph
displaymode = 2;
% vary departure time (false) or vary arrival time (true)
mode = true;


% starting latitude and longitude
phi1 = 0 * (pi / 180);
lambda1 = 0 * (pi / 180);
% target latitude and longitude
phi2 = 0 * (pi / 180);
lambda2 = 180 * (pi / 180);
% conversion to radii
r_prime1 = latlong(phi1, lambda1, a_body, b_body, e_first);
r_prime2 = latlong(phi2, lambda2, a_body, b_body, e_first);


% initial mission time for first possible launch (Gregorian time)
missiont0 = [964 0];
t0 = 0;
% spacing between each examined launch times in seconds
writeinterval = 200;
% examined departure & arrival dates in seconds after t0 (must be a positive multiple of the write interval)
max_dep = 1000;
max_arr = 6000;
% minimum time of flight examined in seconds (must be a positive multiple of the write interval)
tof_min = 200;
% maximum & minimum delta-v plotted
maxdv = 10^6;
mindv = 1;
% maximum & minimum C3 plotted
maxC3 = 100;
minC3 = -100;
ellipsoidres = 2^10;
trajectoryres = 2^14;
bounds = a_body * 2;


if mode == true
    maxcolumn = ((max_arr - tof_min) / (writeinterval)) + 1;
else
    maxcolumn = ((max_dep) / (writeinterval)) + 1;
end
out1 = zeros(1, maxcolumn);
if or(((displaymode) == 2), ((displaymode) == 5))
    out2 = zeros(1, maxcolumn);
end
S_t1 = zeros(6, maxcolumn);
plot_mat = false(1, maxcolumn);
columnindex = 1;
t_dg = missiont0;
[r1, v1, ~] = ECEF2ECI(r_prime1, [0; 0; 0], [0; 0; 0], t_dg, t0, omega_0, sec_matrix, spec_matrix);
S1 = [r1; v1];
tof = tof_min;
if mode == true
while true
    toa = t0 + tof;
    daysa = 0;
    if toa > 86400.002
        Cap = true;
        while Cap == true
            if toa > 86400.002
                toa = toa - 86400.002;
                daysa = daysa + 1;
            else
                Cap = false;
            end
        end
    end
    t_ag = DTG(t_dg, daysa);
    [r2, v2, ~] = ECEF2ECI(r_prime2, [0; 0; 0], [0; 0; 0], t_ag, toa, omega_0, sec_matrix, spec_matrix);
    S2 = [r2; v2];
    if or(((displaymode) == 2), ((displaymode) == 5))
        [plotted, output1, output2, v_out] = Lambert(S1, S2, tof, mu, a_body, e_first, selectionmode, weight, displaymode, epsilon_house1, epsilon_house2, epsilon_halley, epsilon_geoalt, max_house, max_halley, max_geoalt);
        if and((output1 > maxdv), (displaymode == 2)) == true
            output1 = 0;
        end
        if displaymode == 5
            dist = norm(r1);
            energy_escape = (2 * mu) / (dist);
            output1 = ((output1)^2 - energy_escape) * 10^(-6);
            if output1 > maxC3
                output1 = 0;
            end
        end
        if output2 > maxdv
            output2 = 0;
        end
        out1(1, columnindex) = output1;
        out2(1, columnindex) = output2;
        S_t1(:, columnindex) = [r1; v_out];
        if plotted == true
            plot_mat(1, columnindex) = plotted;
        end
    else
        [plotted, output1, ~, v_out] = Lambert(S1, S2, tof, mu, a_body, e_first, selectionmode, weight, displaymode, epsilon_house1, epsilon_house2, epsilon_halley, epsilon_geoalt, max_house, max_halley, max_geoalt);
        if and((output1 > maxdv), ((displaymode == 6) == false)) == true
            output1 = 0;
        end
        if displaymode == 6
            dist = norm(r1);
            energy_escape = (2 * mu) / (dist);
            output1 = ((output1)^2 - energy_escape) * 10^(-6);
            if output1 > maxC3
                output1 = 0;
            end
        end
        out1(1, columnindex) = output1;
        S_t1(:, columnindex) = [r1; v_out];
        if plotted == true
            plot_mat(1, columnindex) = plotted;
        end
    end
    columnindex = columnindex + 1;
    if columnindex > maxcolumn
        break
    end
    tof = tof + writeinterval;
end
end

k = find(plot_mat);
maxcolumn = size(k, 2);
if (maxcolumn > 0) == false
    fprintf('No Valid Trajectories')
end
S_t1old = S_t1;
S_t1 = zeros(6, maxcolumn);
for i = 1:maxcolumn
    search = k(i);
    S_t1(:, i) = S_t1old(:, search);
end
i = 1;
r_mat = zeros(3, (trajectoryres + 1), maxcolumn);
t_mat = tof_min:writeinterval:max_arr;
while true
    S_out = S_t1(:, i);
    Kep_temp = CTK(S_out, mu);
    secondindex = 1;
    delta_t = 0;
    range_t = t_mat(k(1, i));
    while true
        if mode == false
            % assign t0 and t_dg
        end
        S_col = KTC(Kep_temp, mu, delta_t, 0, epsilon_newton, max_newton);
        tc = t0 + delta_t;
        daysa = 0;
        if tc > 86400.002
            Cap = true;
            while Cap == true
                if tc > 86400.002
                    tc = tc - 86400.002;
                    daysa = daysa + 1;
                else
                    Cap = false;
                end
            end
        end
        t_cg = DTG(t_dg, daysa);
        [r_col, ~, ~] = ECI2ECEF(S_col(1:3), [0; 0; 0], [0; 0; 0], t_cg, tc, omega_0, sec_matrix, spec_matrix);
        r_mat(:, secondindex, i) = r_col;
        secondindex = secondindex + 1;
        if secondindex > (trajectoryres + 1)
            break
        end
        delta_t = delta_t + (range_t) / (trajectoryres);
    end
    i = i + 1;
    if i > maxcolumn
        break
    end
end
if and((phi1 == phi2), (lambda1 == lambda2)) == true
    if dot((r_prime1 / norm(r_prime1)), [1; 0; 0]) == 1
        plane_norm = cross(r_prime1, [0; 1; 0]);
    else
        plane_norm = cross([1; 0; 0], r_prime1);
    end
else
    plane_norm = cross(r_prime1, r_prime2);
    plane_norm = (plane_norm) / (norm(plane_norm));
end
y_vec = (r_prime1 / (norm(r_prime1)));
x_vec = cross((y_vec), (plane_norm));
theta_mat = 0:(2 * pi)/(ellipsoidres):(2 * pi);
perimeter = [cos(theta_mat); sin(theta_mat)];
point_mat = perimeter(1) .* x_vec + perimeter(2) .* y_vec;
a_mat = dot(((point_mat) .^(2)), ([a_body; a_body; b_body] .^(-2)));
a_mat = abs(a_mat);
n_mat = sqrt(a_mat) ./ (a_mat);
perimeter = n_mat .* perimeter;
plot(perimeter(1, :), perimeter(2, :), 'blue')
T = [x_vec y_vec];
flat_mat = zeros(2, (trajectoryres + 1), maxcolumn);
plane_norm = repmat(plane_norm, [1 (trajectoryres + 1)]);
for i = 1:maxcolumn
    r_temp = r_mat(:, :, i);
    dist_mat = dot(plane_norm, r_temp) / norm(plane_norm);
    proj_3d = r_temp - ((plane_norm) .* (dist_mat));
    proj_2d = lsqminnorm(T, proj_3d);
    x = proj_2d(1, :);
    y = proj_2d(2, :);
    hold on
    plot(x(:), y(:))
end
r1_dist = dot(plane_norm(:, 1), r_prime1) / norm(plane_norm(:, 1));
r1_3d = r_prime1 - ((plane_norm(:, 1)) .* (r1_dist));
r1_2d = lsqminnorm(T, r1_3d);
r2_dist = dot(plane_norm(:, 1), r_prime2) / norm(plane_norm(:, 1));
r2_3d = r_prime2 - ((plane_norm(:, 1)) .* (r2_dist));
r2_2d = lsqminnorm(T, r2_3d);
plot(r1_2d(1), r1_2d(2), 'o-');
plot(r2_2d(1), r2_2d(2), 'o-');
horizonx = [-bounds; bounds];
horizonslope = cross([r2_2d; 0], [0; 0; 1]);
horizonslope = horizonslope(2) / horizonslope(1);
horizony = zeros(2, 1);
horizony(1) = (-bounds - r2_2d(1)) * horizonslope + r2_2d(2);
horizony(2) = (bounds - r2_2d(1)) * horizonslope + r2_2d(2);
plot(horizonx, horizony, 'red');
axis([-bounds bounds -bounds bounds])
pbaspect([1 1 1])
title('Suborbital Trajectories')
xlabel('metres')
ylabel('metres')
figure;
t_mattemp = t_mat;
t_mat = zeros(1, maxcolumn);
out1temp = out1;
out1 = zeros(1, maxcolumn);
if or((displaymode == 2), (displaymode == 5)) == true
    out2temp = out2;
    out2 = zeros(1, maxcolumn);
end
for i = 1:maxcolumn
    search = k(i);
    t_mat(1, i) = t_mattemp(1, search);
    out1(1, i) = out1temp(1, search);
    if or((displaymode == 2), (displaymode == 5)) == true
        out2(1, i) = out2temp(1, search);
    end
end
if displaymode == 1
    plot(t_mat, out1)
    title('Suborbital Hop Combined Burns')
    ylabel('{\Sigma}{\Delta}v Cost (m/s)')
elseif displaymode == 2
    title('Suborbital Hop Separate Departure & Arrival Burns')
    yyaxis left
    plot(t_mat, out1)
    ylabel('Departure {\Delta}v Cost (m/s)')
    yyaxis right
    plot(t_mat, out2)
    ylabel('Arrival {\Delta}v Cost (m/s)')
elseif displaymode == 3
    plot(t_mat, out1)
    title('Suborbital Hop Departure Burn')
    ylabel('Departure {\Delta}v Cost (m/s)')
elseif displaymode == 4
    plot(t_mat, out1)
    title('Suborbital Hop Arrival Burn')
    ylabel('Arrival {\Delta}v Cost (m/s)')
elseif displaymode == 5
    title('Suborbital Hop Separate Departure & Arrival Burns')
    yyaxis left
    plot(t_mat, out1)
    ylabel('Departure C3 Cost (MJ)')
    yyaxis right
    plot(t_mat, out2)
    ylabel('Arrival {\Delta}v Cost (m/s)')
elseif displaymode == 6
    plot(t_mat, out1)
    title('Suborbital Hop Departure Burn')
    ylabel('Departure C3 Cost (MJ)')
end
xlabel('t+ Arrival Time in Seconds') 




function r_prime = latlong(phi, lambda, a_body, b_body, e_first)
if phi == pi / 2 
    z = b_body;
elseif phi == -pi / 2
    z = -b_body;
elseif phi == 0
    z = 0;
else
    z = a_body * (1 - e_first^2) * sin(phi) / (sqrt(1 - e_first^2 * (sin(phi))^2));
end
if abs(phi) == pi / 2
    p = 0;
elseif phi == 0
    p = a_body;
else
    p = a_body * cos(phi) / (sqrt(1 - e_first^2 * (sin(phi))^2));
end
x = p * cos(lambda);
y = p * sin(lambda);
r_prime = [x; y; z];
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


function [plotted, output1, output2, v_out] = Lambert(S1, S2, t, mu, a_body, e_first, selectionmode, weight, displaymode, epsilon_house1, epsilon_house2, epsilon_halley, epsilon_geoalt, max_house, max_halley, max_geoalt)
arguments
    S1
    S2
    t {mustBeGreaterThan(t, 0)}
    mu {mustBeGreaterThan(mu, 0)}
    a_body
    e_first
    selectionmode {mustBeInteger, mustBeInRange(selectionmode, 1, 4)}
    weight
    displaymode {mustBeInteger, mustBeInRange(displaymode, 1, 6)}
    epsilon_house1
    epsilon_house2
    epsilon_halley
    epsilon_geoalt
    max_house {mustBeInteger}
    max_halley {mustBeInteger}
    max_geoalt {mustBeInteger}
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
VOut = zeros(3, matrixsize);
VIn = zeros(3, matrixsize);
e_mat = zeros(1, matrixsize);
for i = 1:matrixsize
    x = SolutionsList(i, 1);
    y = SolutionsList(i, 2);
    vr1 = gamma * (((lambda * y) - x) - (rho * ((lambda * y) + x))) / r1norm;
    vr2 = -gamma * (((lambda * y) - x) + (rho * ((lambda * y) + x))) / r2norm;
    vt1 = gamma * sigma * (y + (lambda * x)) / r1norm;
    vt2 = gamma * sigma * (y + (lambda * x)) / r2norm;
    vdep = vr1 * a1 + vt1 * b1;
    varr = vr2 * a2 + vt2 * b2;
    VOut(:, i) = vdep;
    VIn(:, i) = varr;
    Kepler_temp = CTK([r1; vdep], mu);
    e_mat(1, i) = Kepler_temp(2);
end
ROut = repmat(r1, [1 matrixsize]);
v1 = repmat(v1, [1 matrixsize]);
geo_norm = zeros(3, matrixsize);
for i = 1:matrixsize
geo_norm(:, i) = geonorm(ROut(:, i), a_body, e_first, epsilon_geoalt, max_geoalt);
end
D = dot(geo_norm, (VOut - v1));
D = (D >= 0);
e_mat = (e_mat < 1);
D = ((D) & (e_mat));
k = find(D);
matrixsize = size(k, 2);
if matrixsize > 0
    plotted = true;
    v1 = repmat(v1, 1, matrixsize);
    v2 = repmat(v2, 1, matrixsize);
    VelocityDeparture = zeros(3, matrixsize);
    VelocityArrival = zeros(3, matrixsize);
    for i = 1:matrixsize
        search = k(i);
        VelocityDeparture(:, i) = VOut(:, search);
        VelocityArrival(:, i) = VOut(:, search);
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
    elseif displaymode == 2
        output1 = norm(RelVelDep(:, Index));
        output2 = norm(RelVelArr(:, Index));
    elseif displaymode == 3
        output1 = norm(RelVelDep(:, Index));
        output2 = 0;
    elseif displaymode == 4
        output1 = norm(RelVelArr(:, Index));
        output2 = 0;
    elseif displaymode == 5
        output1 = norm(VelocityDeparture(:, Index));
        output2 = norm(RelVelArr(:, Index));
    elseif displaymode == 6
        output1 = norm(VelocityDeparture(:, Index));
        output2 = 0;
    end
    v_out = VelocityDeparture(:, Index);
else
    plotted = false;
    output1 = 0;
    output2 = 0;
    v_out = [0; 0; 0];
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


function [r, v, a] = ECEF2ECI(r_prime, v_prime, a_prime, GreTime, t, omega_0, sec_matrix, spec_matrix)
[rot, rot_dot, rot_doubledot] = rotfunction(GreTime, t, omega_0, sec_matrix, spec_matrix);
T = anglefunction(rot);
omega = omegafunction(rot, rot_dot);
omega_dot = omegadotfunction(rot, rot_dot, rot_doubledot);
r = T * (r_prime);
v = T * (v_prime) + cross((omega), (r_prime));
a = T * (a_prime) + (cross((omega_dot), (r_prime))) + (2 * (cross((omega), (v_prime)))) + (cross((omega), (cross((omega), (r_prime)))));
end


function [r_prime, v_prime, a_prime] = ECI2ECEF(r, v, a, GreTime, t, omega_0, sec_matrix, spec_matrix)
[rot, rot_dot, rot_doubledot] = rotfunction(GreTime, t, omega_0, sec_matrix, spec_matrix);
T = anglefunction(rot);
omega = omegafunction(rot, rot_dot);
omega_dot = omegadotfunction(rot, rot_dot, rot_doubledot);
r_prime = lsqminnorm(T, r);
v_prime = lsqminnorm(T, v) - cross((omega), (r_prime));
a_prime = lsqminnorm(T, a) - (cross((omega_dot), (r_prime))) - (2 * (cross((omega), (v_prime)))) - (cross((omega), (cross((omega), (r_prime)))));
end

% transformation matrix to convert radii in the ECEF frame to radii in the ECI frame
function rottransform = anglefunction(rot)
psi = rot(1);
theta = rot(2);
omega = rot(3);
precession = [1 0 0; 0 cos(psi) -sin(psi); 0 sin(psi) cos(psi)];
nutation = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
rotation = [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];
rottransform = (rotation) * ((nutation) * (precession));
end


% angular velocity of the ECEF frame in the ECI frame
function omega = omegafunction(rot, rot_dot)
psi = rot(1);
theta = rot(2);
psi_dot = rot_dot(1);
theta_dot = rot_dot(2);
omega_dot = rot_dot(3);
omega_x = ((omega_dot) * sin(theta) * sin(psi)) + ((theta_dot) * cos(psi));
omega_y = ((-(omega_dot)) * sin(theta) * cos(psi)) + ((theta_dot) * sin(psi));
omega_z = ((omega_dot) * cos(theta)) + (psi_dot);
omega = [omega_x; omega_y; omega_z];
end


% angular acceleration of the ECEF fraction in the ECI frame
function omega_dot = omegadotfunction(rot, rot_dot, rot_doubledot)
psi = rot(1);
theta = rot(2);
psi_dot = rot_dot(1);
theta_dot = rot_dot(2);
omega_dot = rot_dot(3);
psi_doubledot = rot_doubledot(1);
theta_doubledot = rot_doubledot(2);
omegadot_x = ((theta_dot) * ((omega_dot) * cos(theta) * sin(psi))) + ((theta_doubledot) * cos(psi)) + ((psi_dot) * (((omega_dot) * sin(theta) * cos(psi)) - ((theta_dot) * sin(psi))));
omegadot_y = ((theta_dot) * ((-(omega_dot)) * cos(theta) * cos(psi))) + ((theta_doubledot) * sin(psi)) + ((psi_dot) * (((omega_dot) * sin(theta) * sin(psi)) + ((theta_dot) * cos(psi))));
omegadot_z = ((theta_dot) * ((-(omega_dot)) * sin(theta))) + (psi_doubledot);
omega_dot = [omegadot_x; omegadot_y; omegadot_z];
end


function [rot, rot_dot, rot_doubledot] = rotfunction(GreTime, t, omega_0, sec_matrix, spec_matrix)
psi_vec      = sec_matrix(1, :);
theta_vec    = sec_matrix(2, :);
psi_ipvec    = spec_matrix(:, 1);
psi_oopvec   = spec_matrix(:, 2);
theta_ipvec  = spec_matrix(:, 3);
theta_oopvec = spec_matrix(:, 4);
omega_vec    = spec_matrix(:, 5);
chi_vec      = spec_matrix(:, 6);
Plus900 = GTE(GreTime) + t;
t_reduced = Plus900 / ((10^4) * 365.25 * 86400);
secpol = [1; t_reduced; (t_reduced)^2; (t_reduced)^3; (t_reduced)^4; (t_reduced)^5; (t_reduced)^6; (t_reduced)^7; (t_reduced)^8; (t_reduced)^9; (t_reduced)^10];
 % accumulated secular precession
psi_a = psi_vec * secpol;
% accumulated secular obliquity
theta_a = theta_vec * secpol;
% argument of spectral nutation
xi_vec = ((Plus900) * (omega_vec)) + (chi_vec);
% spectral nutation
delta_psi = sum(times((psi_ipvec), (sin(xi_vec))) + times((psi_oopvec), (cos(xi_vec))));
delta_theta = sum(times((theta_ipvec), (cos(xi_vec))) + times((theta_oopvec), (sin(xi_vec))));
% time for purposes of determining current point in the sidereal day
SidPhase = SPh(GreTime, t);
% final angles
psi = psi_a + delta_psi;
theta = theta_a + delta_theta;
omega = SidPhase * omega_0;
rot = [psi; theta; omega];


% derivatives
secpol_prime = [0; 1; 2*(t_reduced); 3*(t_reduced)^2; 4*(t_reduced)^3; 5*(t_reduced)^4; 6*(t_reduced)^5; 7*(t_reduced)^6; 8*(t_reduced)^7; 9*(t_reduced)^8; 10*(t_reduced)^9];
% secular derivatives
psi_aprime = (psi_vec * secpol_prime) / ((10^4) * 365.25 * 86400);
theta_aprime = (theta_vec * secpol_prime) / ((10^4) * 365.25 * 86400);
% spectral nutation derivatives
delta_psiprime = sum((times((times((psi_ipvec), (omega_vec))), (cos(xi_vec)))) - (times((times((psi_oopvec), (omega_vec))), (sin(xi_vec)))));
delta_thetaprime = sum((-times((times((theta_ipvec), (omega_vec))), (sin(xi_vec)))) + (times((times((theta_oopvec), (omega_vec))), (cos(xi_vec)))));
% final derivatives
psi_dot = psi_aprime + delta_psiprime;
theta_dot = theta_aprime + delta_thetaprime;
omega_dot = omega_0;
rot_dot = [psi_dot; theta_dot; omega_dot];


% second derivatives
secpol_doubleprime = [0; 0; 2; 6*(t_reduced); 12*(t_reduced)^2; 20*(t_reduced)^3; 30*(t_reduced)^4; 42*(t_reduced)^5; 56*(t_reduced)^6; 72*(t_reduced)^7; 90*(t_reduced)^8];
% secular second derivatives
psi_adoubleprime = (psi_vec * secpol_doubleprime) / (((10^4) * 365.25 * 86400)^2);
theta_adoubleprime = (theta_vec * secpol_doubleprime) / (((10^4) * 365.25 * 86400)^2);
% spectral nutation second derivatives
delta_psidoubleprime = sum((-times((times((psi_ipvec), (power((omega_vec), 2)))), (sin(xi_vec)))) + (-times((times((psi_oopvec), (power((omega_vec), 2)))), (cos(xi_vec)))));
delta_thetadoubleprime = sum((-times((times((theta_ipvec), (power((omega_vec), 2)))), (cos(xi_vec)))) + (-times((times((theta_oopvec), (power((omega_vec), 2)))), (sin(xi_vec)))));
psi_doubledot = psi_adoubleprime + delta_psidoubleprime;
theta_doubledot = theta_adoubleprime + delta_thetadoubleprime;
rot_doubledot = [psi_doubledot; theta_doubledot];
end


function Plus900 = GTE(GreTime)
Yrs = GreTime(1) - 900;
LYrs = floor(Yrs/4);
Yrs = Yrs - ((LYrs) * 4);
Plus900 = ((LYrs) * 4 * 365.25 * 86400.002) + ((Yrs) * 365 * 86400.002) + ((GreTime(2)) * 86400.002);
end


function SidPhase = SPh(GreTime, t)
SidPhase = ((GreTime(2)) * 86400.002) + t;
end


function geo_norm = geonorm(r, a_body, e_first, epsilon_geoalt, max_geoalt)
    x = r(1);
    y = r(2);
    p = [x; y];
    norm_p = sqrt(x^2 + y^2);
    p = (p) / (norm_p);
    z = abs(r(3));
    if norm_p == 0
        if z == 0
            geo_norm = [0; 0; 0];
        else
            if z > 0
                phi = pi;
                p = cos(phi) * p;
                z_norm = sin(phi);
                geo_norm = [p; z_norm];
            else
                phi = -pi;
                p = cos(phi) * p;
                z_norm = sin(phi);
                geo_norm = [p; z_norm];
            end
        end
    else
        c = a_body * (e_first)^2;
        eprime = sqrt(1 - (e_first)^2);
        zprime = eprime * z;
        u = 2 * (zprime - c);
        v = 2 * (zprime + c);
        tM = (c - zprime) / norm_p;
        fM = norm_p * (tM)^4 + u * (tM)^3 + v * (tM) - norm_p;
        t1 = (norm_p - c + zprime) / (norm_p - c + 2 * zprime);
        t0 = norm_p / (zprime + c);
        if tM <= 0
            t = t1;
        elseif tM >= 1 
            t = t0;
        elseif fM >= 0
            t = t0;
        else
            t = t1;
        end
        for i = 1:max_geoalt
            deltat = (norm_p - (norm_p * (t)^4 + u * (t)^3 + v * t)) / (4 * norm_p * (t)^3 + 3 * u * (t)^2 + v);
            lastvalue = t;
            t = t + deltat;
            TE = abs(t - lastvalue);
            if TE <= epsilon_geoalt
                break
            end
        end
        if r(3) > 0
            phi = atan2((1 - t^2), (2 * eprime * t));
        else
            phi = -atan2((1 - t^2), (2 * eprime * t));
        end
        p = cos(phi) * p;
        z_norm = sin(phi);
        geo_norm = [p; z_norm];
    end
end

