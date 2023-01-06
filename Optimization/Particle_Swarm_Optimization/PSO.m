clc;
clear;
close all;
%% Problem Definition
% initCostfunction = @(x) cost_func(x);   % cost function
% penalityFunction = @(x) penalty_func(x);% penality Function
% Cp =@(t) sqrt(t+1);
% Costfunction =@(x,t) cost_func(x) ...   % Objective function to  minimize
%                     + Cp(t)*penalty_func(x);
% nVar =4;                            % Lower Bound of decision variables 
% varSize = [1 nVar];                 % Matix Size of decision Variables
% varMin = [0.125,0.1,0.1,0.125];     % Lower bound of Decision varibles
% varMax = [5,10,10,5];               % Upper bound of devision 
% t = 0;
%% Test function
Costfunction =@(x) Test_func(x);
nVar =2;                            % Lower Bound of decision variables 
varSize = [1 nVar];                 % Matix Size of decision Variables
varMin = [-10 -10];                 % Lower bound of Decision varibles
varMax = [10,10];               	% Upper bound of devision 
t = 0;
%% Parameters of PSO
maxIterations = 100;                % Maximum Number of iteratios
nPopulation = 50;                   % Population Size
wMax = 0.9;
wMin = 0.2;
w = wMax-((wMax-wMin)/100)*t;       % Intertia Coefficient
c1 = 2;                             % cognitive Acceleration Coefficient
c2 = 2.2;                             % social Acceleration Coefficient
% Calulate the constriction factor
Phi = c1 + c2;
if Phi > 4
    kappa = 2/abs(2-Phi-sqrt(Phi^2-4*Phi));
else
    kappa = 1;
end
%% Initialization
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = 0;
empty_particle.Best.Position = [];
empty_particle.Best.Cost = 0;

% Create population struct array
particle = repmat(empty_particle,nPopulation,1);

% Initialize Global Best
GlobalBest.Cost = inf;
% Initialize Population to random uniform values
for i=1:nPopulation
    % Generate Random solutions
    for j = 1:nVar
        particle(i).Position(j) = unifrnd(varMin(j),varMax(j));
    end
    % Intialize
    particle(i).Velocity = zeros(varSize);
    % Cost Evaluation 
    %particle(i).Cost = Costfunction(particle(i).Position,t);
    particle(i).Cost = Costfunction(particle(i).Position);
    % Update Personal best
    particle(i).Best(1).Position = particle(i).Position;
    particle(i).Best(1).Cost = particle(i).Cost;
    % Update Global Best
    if particle(i).Best(1).('Cost') < GlobalBest.('Cost')
        GlobalBest = particle(i).Best;
        %disp("Update Global best in init");
    end
end
t = t + 1;
% Array to hold best cost value from each population
bestCosts = zeros(maxIterations,3);
%% Main loop of PSO
for outer=1:maxIterations
   for i =1:nPopulation
       % Update Velocity
        particle(i).Velocity = kappa*(w*particle(i).Velocity...
            +c1*rand().*(particle(i).Best.Position - particle(i).Position)...
            +c2*rand().*(GlobalBest.Position - particle(i).Position));
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        % Test boundary conditions 
        for j=1:nVar
            if particle(i).Position(j) < varMin(j)
                particle(i).Position(j) = unifrnd(varMin(j),varMax(j));
            elseif particle(i).Position(j) > varMax(j)
                particle(i).Position(j) = unifrnd(varMin(j),varMax(j));
            end
        end
        % Cost Evaulation
        %particle(i).Cost = Costfunction(particle(i).Position,t);
        particle(i).Cost = Costfunction(particle(i).Position);
        % Update Personal best
        if  particle(i).('Cost') < particle(i).Best(1).('Cost')
             particle(i).Best(1).Position = particle(i).Position;
             particle(i).Best(1).('Cost') = particle(i).('Cost');
             %disp("Update personal best in main loop");
             % Update Global Best
             if particle(i).Best(1).Cost < GlobalBest.Cost
                  GlobalBest = particle(i).('Best');
                  %disp("Update Global best in main loop");
             end    
        end
   end
   t = t + 1;
   % store the best Cost value
   bestCosts(outer,1) = GlobalBest.('Cost');
   bestCosts(outer,2) = GlobalBest.Position(1);
   bestCosts(outer,3) = GlobalBest.Position(2);
   % display Iteration info
   disp(['Iteration #:' num2str(outer) ' Best Cost =' num2str(GlobalBest.Cost) ' Best Position =' num2str(GlobalBest.Position)]);
end
%% Results

plot(bestCosts(:,1))
title("Best Cost Vs interation")
subtitle("Booth Test function")
xlabel("Number of Iterations")
ylabel("Best Cost")
legend("Best Cost")
dim = [.67 .7 .1 .1];
str = {"PSO parameters:"," c1="+num2str(c1), " c2="+num2str(c2), " kappa="+num2str(kappa), " w="+num2str(w)};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
 disp(['Best Cost =' num2str(bestCosts(100))]);