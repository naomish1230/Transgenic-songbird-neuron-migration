close all; clear all;

%% Setting up 3D visual space for simulation

% Creating 3D sphere with a 1 mm diameter to approximate the shape of HVC
[X,Y,Z] = sphere; 
X = X+1; Y = Y+1; Z = Z+1; X = X*0.5; Y = Y*0.5; Z = Z*0.5; 
surf(X,Y,Z,'FaceColor', 'interp','FaceAlpha',0.2,'EdgeColor','none'); 
axis square; axis equal; grid off;

hold on

% Creating semi-transparent ventricular zone (VZ) plane to intersect the spherical HVC
VZplaneX = [0 0 1 1]; VZplaneY = [0 1 1 0]; VZplaneZ = [0.51 0.51 0.85 0.85]; 
patch(VZplaneX, VZplaneY, VZplaneZ,'black'); alpha 0.05;


% For aesthetic purposes, we trim this spherical HVC to not extend beyond of the LV plane
for trimcount = 1:length(X)
    for trimcount2 = 1:length(X)
        if Z(trimcount,trimcount2) > 0.51 + 0.34*(X(trimcount,trimcount2))
        Z(trimcount,trimcount2)= 0.51 + 0.34*(X(trimcount,trimcount2));
        end
    end
end

%% Defining in vivo-derived distributions for generating movement

LND_stepdist = makedist('Lognormal','mu',0.83367,'sigma',0.79433); % Step size (d) distribution
logdist_thresholded_thresh5 =makedist('Logistic','mu', -1.37616,'sigma',9.79152); % Angle (a) change distribution

%% User-Defined Parameters for Simulation Run

% User decides how long simulation is run for and how many cells are being simulated
num_steps = 1008; % Number of time steps that simulation is run for, with each step being equivalent to 10 real world minutes. For example, 7 days is 1008 steps.
num_cells = 100;

rsquare_thresh = 0.95; % Adjusted R^2 threshold for determing goodness of fit of curves fit to Mean Squared Displacement (MSD) that the alpha parameter is derived from
plot_traces = 1; % Determines whether traces of migration are shown: 0 is off, 1 is on

%% Simulating Birth and Migration of Cells

for cell_ID = 1:num_cells
    
    % First, cell must be born at some random point in 3D within HVC on the VZ plane facing in some random direction
    [x(1),y(1),z(1)] = birth_cell(); 
    plot3(x(1),y(1), z(1),'g.','MarkerSize',10) % Plotting birth point

    starttheta(cell_ID,:) = 360*rand(); % Initial angle that cell is facing in XY
    startphi(cell_ID,:) = 360*rand(); % Initial angle that cell is facing in XZ
        
    for steps = 1:num_steps % Generating movement in a stepwise manner        
        
        dist_prob = rand(); % Creating random probability to determine which angle distribution is drawn from for the following step

        if dist_prob <= 0.9 
            theta(cell_ID,steps) = random(logdist_thresholded_thresh5);
            phi(cell_ID,steps) = random(logdist_thresholded_thresh5);
        else
            theta(cell_ID,steps) =unifrnd(-180,180);
            phi(cell_ID,steps) = unifrnd(-180,180);
        end

        if steps == 1
            nexttheta(cell_ID,steps) = starttheta(cell_ID,:) + theta(cell_ID,steps);
            nextphi(cell_ID,steps) = startphi(cell_ID,:) + phi(cell_ID,steps);
        else
            nexttheta(cell_ID,steps) = theta(cell_ID,steps) + nexttheta(cell_ID,steps-1);
            nextphi(cell_ID,steps) = phi(cell_ID,steps) + nextphi(cell_ID,steps-1);
        end
        
        step_dist(cell_ID,steps) = random(LND_stepdist)/1000; % Drawing random step distance from distribution from in vivo data and converting from microns to millimeters
       
        dx = step_dist(cell_ID,steps)*cosd(nextphi(cell_ID,steps))*sind(nexttheta(cell_ID,steps));
        dy = step_dist(cell_ID,steps)*cosd(nextphi(cell_ID,steps))*cosd(nexttheta(cell_ID,steps));
        dz = step_dist(cell_ID,steps)*sind(nextphi(cell_ID,steps));

        x(steps + 1) = x(steps) + dx;
        y(steps + 1) = y(steps) + dy;
        z(steps + 1) = z(steps) + dz;

        % Encoding for condition that cells do not migrate back into VZ post-birth. If cell attempts to do so, their angle and step size are re-drawn from in vivo-derived distributions
        while z(steps + 1) >= 0.51 + (0.34*x(steps+1)) 
            if dist_prob <= 0.9
                    theta(cell_ID,steps) =random(logdist_thresholded_thresh5);
                    phi(cell_ID,steps) = random(logdist_thresholded_thresh5);
                    else
                    theta(cell_ID,steps) =unifrnd(-180,180);
                    phi(cell_ID,steps) = unifrnd(-180,180);
                    end
            
            if steps == 1
                    nexttheta(cell_ID,steps) = starttheta(cell_ID,:) + theta(cell_ID,steps);
                    nextphi(cell_ID,steps) = startphi(cell_ID,:) + phi(cell_ID,steps);
            else
                    nexttheta(cell_ID,steps) = theta(cell_ID,steps) + nexttheta(cell_ID,steps-1);
                    nextphi(cell_ID,steps) = phi(cell_ID,steps) + nextphi(cell_ID,steps-1);
            end
          
            step_dist(cell_ID,steps) = random(LND_stepdist)/1000; % Drawing random step distance from distribution from in vivo data and converting from microns to millimeters

            dx = step_dist(cell_ID,steps)*cosd(nextphi(cell_ID,steps))*sind(nexttheta(cell_ID,steps));
            dy = step_dist(cell_ID,steps)*cosd(nextphi(cell_ID,steps))*cosd(nexttheta(cell_ID,steps));
            dz = step_dist(cell_ID,steps)*sind(nextphi(cell_ID,steps));

            x(steps + 1) = x(steps) + dx;
            y(steps + 1) = y(steps) + dy;
            z(steps + 1) = z(steps) + dz;
        end

        store_dist(cell_ID,steps) = step_dist(cell_ID,steps); % Storing distance traveled in a step for later, each row is a cell each column is a time step
    
    end 

    % At this point, one cell has completed its migration and can have its speed, tortuosity, and MSD alpha parameter calculated
    plot3(x(end),y(end),z(end),'r.','MarkerSize',10) % Plotting endpoint of this cell's migration
    
    % Plotting trace of this cell's migration (when selected by user)
    if plot_traces == 1
        plot3(x,y,z,'k-') 
    end

    % Calculating average speed of cell's migration in um/hr 
    avgspeed(cell_ID,:) = (sum(store_dist(cell_ID,:),2)/((num_steps*10))*(60*1000));
    
    % Calculating tortuosity of cell's migration. If simulation is run for more than 3 hours, a sliding 3 hour average window is used to calculate average tortuosity across entirety of migration
    if num_steps > 18 
        for slideave = 1:(num_steps-17) 
            sliding_tort(slideave,1) = (sum(store_dist(slideave:slideave +17),2))/sqrt(((x(slideave+17)-x(slideave))^2) + ((y(slideave+17) - y(slideave))^2) + ((z(slideave+17) - z(slideave))^2));
        end
        tort(cell_ID,:) = mean(sliding_tort);
    else
        tort(cell_ID,:) = (sum(store_dist(cell_ID,:),2))/sqrt(((x(end)-x(1))^2) + ((y(end) - y(1))^2) + ((z(end) - z(1))^2)); %tortuosity of each cell's migration, total distance traveled/total displacement
    end
        

   % Calculating 2D MSD and deriving alpha parameter. 2D MSD is used over 3D MSD to remain consistent with in vivo tracking experiments but 2D and 3D MSD are generally similar.

   timelag = [1:num_steps]'; % Creating horizontal vector for use in fitting curve to MSD during migration
    
   [MSD_3D, MSD_2D] = MSD_func(x,y,z); 
   
   % Deriving alpha parameter from 2D MSD
   [f1,MSD2Dgof] = fit(timelag,MSD_2D,'power1'); %2D MSD
   Allf1_2D{cell_ID,:} = f1;
   AllMSD2Dgof{cell_ID,:} = MSD2Dgof;
   AllRsquare2D(cell_ID,1) = MSD2Dgof.adjrsquare;
   AllMSD_2D{cell_ID,:} = MSD_2D; %storing MSD from each cell
    
   if MSD2Dgof.adjrsquare >= rsquare_thresh % This limits alpha params. used in analysis to those that originate from curves that are user-defined as good fits as determined by R^2
        Allalpha_2D(cell_ID,1) =f1.b;
   else
        Allalpha_2D(cell_ID,1) = NaN;
   end
    
   % Deriving alpha parameter from 3D MSD
   [f1,MSD3Dgof] = fit(timelag,MSD_3D,'power1');
   Allf1_3D{cell_ID,:} = f1;
   AllMSD3Dgof{cell_ID,:} = MSD3Dgof;
   AllRsquare3D(cell_ID,1) = MSD3Dgof.adjrsquare;
   AllMSD_3D{cell_ID,:} = MSD_3D; 
    
   if MSD3Dgof.adjrsquare >= rsquare_thresh % This limits alpha params. used in analysis to those that originate from curves that are user-defined as good fits as determined by R^2
       Allalpha_3D(cell_ID,1) =f1.b;
   else
       Allalpha_3D(cell_ID,1) = NaN;
   end
end

% Calculating population MSDs for entire cohort of simulaed cells

% 2D population MSD
everyMSD2D = cell2mat(AllMSD_2D);
reshapeMSD2D = reshape(everyMSD2D,[length(MSD_2D),num_cells]);
aveMSD2D = mean(reshapeMSD2D,2);
[f_aveMSD,aveMSD_gof] = fit(timelag,aveMSD2D,'power1');
aveAlpha2D = f_aveMSD.b;

% 3D population MSD
everyMSD3D = cell2mat(AllMSD_3D);
reshapeMSD3D = reshape(everyMSD3D,[length(MSD_3D),num_cells]);
aveMSD3D = mean(reshapeMSD3D,2);
[f_aveMSD,aveMSD_gof] = fit(timelag,aveMSD3D,'power1');
aveAlpha3D = f_aveMSD.b;

%% Plotting Outcome of Simulated Migration

xlabel('X (mm)'); ylabel('Y (mm)');zlabel('Z (mm)');

days = round((num_steps*10)/(60*24)); % Converting length of simulation into days
if days >= 2
    title(sprintf('3D Simulation of %.0d Cells in HVC over %.0d Days',num_cells,days))
else
    hours = steps*10/60;
    title(sprintf('3D Simulation of %.0d Cells in HVC over %0.2f Hours',num_cells,hours))
end

%% Functions 
function [x1,y1,z1] = birth_cell() 
% This function will randomly generate the birth point of a cell.
% First, a random point in the X dimension that is on the VZ plane intersecting with HVC is generated
    x(1) = rand();
    while x(1) > 0.8902 % Checking that randomly generated X position is not outside of HVC to begin with. If outside of HVC, X position is re-generated
         x(1) = rand();
    end
% Next, the range of acceptable Y positions based on this X coordinate is calculated. A random value from this range is then selected to be the Y birth position. 
    yint1 = abs((sqrt(0.25-(x(1)-0.5)^2))- 0.5);
    yint2 = (sqrt(0.25-(x(1)-0.5)^2))+0.5; 
    y(1) = unifrnd(yint1,yint2); 
% Finally, the range of acceptable Z positions based on the previously generated X and Y coordinates is calculated. 

    zint1 = abs((sqrt(0.25-((x(1)-0.5)^2)-((y(1)-0.5)^2)))- 0.5); 
    zint2 = (sqrt(0.25-((x(1)-0.5)^2)-((y(1)-0.5)^2)))+ 0.5; 
 % The corresponding Z position is calculated based on the previously generated X coordinate. If this Z value is outside of the previously calculated range, X and Y are re-generated 
    z(1) = 0.51 + (0.34*x(1)); 
    while z(1) < zint1 || z(1) > zint2 
    x(1) = rand(); 
        while x(1) > 0.8902 
            x(1) = rand();
        end

        yint1 = abs((sqrt(0.25-(x(1)-0.5)^2))- 0.5); 
        yint2 = (sqrt(0.25-(x(1)-0.5)^2))+0.5; 
        y(1) = unifrnd(yint1,yint2);
        z(1) = 0.51 + (0.34*x(1)); 
        zint1 = abs((sqrt(0.25-((x(1)-0.5)^2)-((y(1)-0.5)^2)))- 0.5); 
        zint2 = (sqrt(0.25-((x(1)-0.5)^2)-((y(1)-0.5)^2)))+ 0.5; 
    end
    x1 = x(1); y1 = y(1); z1 = z(1);
end
    
function [MSD_3D,MSD_2D] = MSD_func(x,y,z)
% This function will calculate the individual 2D and 3D MSD for each cell's migration
for tau = 1:length(x)-1
    MSDdx = x(1:1:(end-tau)) - x((tau+1):1:end); 
    MSDdy = y(1:1:(end-tau)) - y((tau+1):1:end);
    MSDdz = z(1:1:(end-tau)) - z((tau+1):1:end);
    
    dr2_3D = MSDdx.^2 +MSDdy.^2 + MSDdz.^2; 
    dr2_2D = MSDdx.^2 +MSDdy.^2; 

    MSD_3D(tau,1) = mean(dr2_3D); % Calculating the mean of the squared displacement in 3D
    MSD_2D(tau,1) = mean(dr2_2D); % Calculating the mean of the squared displacement in 2D
end
end




