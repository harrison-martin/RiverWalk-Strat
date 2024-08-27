%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caitlin Sifuentes - River Avulsion Stratigraphy                         %
% csifuen@iu.edu                                                          %
% March  2021                                                             %
% Goal: Add stratigraphy (3D of RiverWalk) and analyze the spatial        %
% distributions of the channels in cross section                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reset
clearvars 
close all
%% Timer
tic
%% Input Variables
numBlocks = 300; % how many mat files/blocks do we have? ---> might clean this up using names=dir('*.mat');
blockSize = 1000; %number of time steps in each stratigraphic block 
numSteps = numBlocks*blockSize; % total thickness/number of timesteps
dt = 10; % timestep length (yrs)
tFinal = numSteps*dt; % total runtime (yrs)
dx=500; %cell size in meters
domainWidth = 301; % number of cells (cell = 500m)
stratstart = 5; %what timestep to start on...

domainDistance = 301; % used to initiate timestep and facies strat
res = 10; % cell/m, res = 10 is a cell resolution 10 cm 
res_m = 1/res; %mcell thickness in meters
totalRunNum = 3; % or 5, used for creating modal percentages

% batch run inputs
toggle_input = [0, 0, 0, 0, 0, 0, 0, 0, 0]; % choose what type of cross section you want to take! 0 = strike/width, 1 = dip/length
sliceLocal_input = [25, 50, 75, 270, 260, 250, 200, 175, 150]; % along which row or column do you want to take a slice of the basin? 
ChiPos = {'Chi25','Chi50','Chi75','270','260','250','200','175','150'}; %% BACKWARDS []
% TitlingForRuns2 = {'/RW_ApexAvul','/SD_ApexAvul','/RW_MultiAvulLocal','/SD_MultiAvulLocal','/AttractionOnly','/RepulsionOnly','/AttractandRepulsTEST','/DepoHeal','/FarfieldHeal'};
TitlingForRuns = {'/RW_ApexAvul','/SD_ApexAvul','/RW_MultiAvulLocal','/SD_MultiAvulLocal','/AttractionOnly','/RepulsionOnly','/AttractandRepuls','/DepoHeal','/FarfieldHeal'};
HealDirectionForRuns = [2,2,2,2,2,2,2,2,2]; % choice in healing direction: 0 deposition-only, 1 erosion-only, 2 mid-directed, 3 far-field directed
FASforRuns = [3,3,3,3,3,3,3,3,3];

load faciesColors

for runNumTotal = 3
        exportDirectory = strcat('/Volumes/G1/CaitlinSifuentes_MSFiles/HCH2/FAS3/',TitlingForRuns{runNumTotal},'/Data');
        exportDirectory2 = strcat('/Volumes/G1/CaitlinSifuentes_MSFiles/HCH2/FAS3/',TitlingForRuns{runNumTotal},'/Figures');
        cd(exportDirectory)
        
    for runNum = 2 % make sure the var totalRunNum is the same as the value for runNum
        % load batch run vars
        toggle = toggle_input(runNum);
        sliceStruct = load(['Chi',num2str(sliceLocal_input(runNum)),'.mat']); 
        sliceLocal = getfield(sliceStruct,strcat('hardcodeChi',num2str(sliceLocal_input(runNum)),'Local'));
        exportName = ChiPos{runNum}; % used to save figures

    %% Load Topographic Elevations and Elevation Differences
    loadedStrike = load(['StrikeTracker',num2str(runNum),'.mat']);
    load('StartingSurface.mat')
    
    if runNum == 1
        thisStrike = loadedStrike.strikeTracker1;
    elseif runNum == 2
        thisStrike = loadedStrike.strikeTracker2;
    elseif runNum == 3
        thisStrike = loadedStrike.strikeTracker3;
    else
        error('Invalid runNum... only works with 3 Chi''s right now')
    end

    %undo subsidence
    load('Subsidence.mat')
    tsu = 1:numSteps;
    tsu=tsu';
    thisStrike(tsu,:) = thisStrike(tsu,:) + dt*tsu*betaSubMtx(sliceLocal,1); %move the surface back up to correct for subsidence 
    
    % Create Basin Volume Slice (2D)
    elevDiff = zeros(size(thisStrike));
    elevDiff(2:numSteps,:) = thisStrike(2:numSteps,:) - thisStrike(1:numSteps-1,:); % Use for loop to look at differences between layers to see if that^^^works, should be the exact same as totalElevDiffSlice
    
    %% Facies Matrix (facies)
    % Initiate totalFaciesSlice
    %totalFaciesSlice = zeros(numSteps,domainDistance);

        for i = 1:numBlocks
            i
           % load in one block at a time
            load(['FaciesMtx',num2str(i),'.mat']); % loads current block "i" as a double
    
            % take a current slice of the block and add it to the total slice
            if toggle == 0
                wxSection3DFacies = facies(sliceLocal,:,:); % creates a cross section at the chosen row (width parallel x-section)
                currentFaciesSlice = squeeze(wxSection3DFacies)'; % returns an array with the same elements as the input array, but with dimensions of length 1 removed resulting in a 2D matrix
                m = (1:blockSize:numSteps); % used to index where new slices of currentElevDiffSlice are placed
                totalFaciesSlice(m(i):blockSize*i,1:domainWidth) = currentFaciesSlice; % fills totalElevDiffSlice with the currentElevDiffSlice
            end
    
            % clear 3D blocks
            clear facies   
        end
    %load('Strike151_totalFacies') %this was made running the commented code above, but that code takes a long time, so save and dont keep running it

%% Stratigraphy using Volume(timestepStrat & faciesStrat)
    % Creating Stratigraphy in Meters 

    % FACIES & TIMESTEP MATRIX
    % Initiation of Stratigraphic Cross Section
     if toggle == 0
         stratCrossSec = thisStrike; % use 2D newSurf for stratigraphic cross section
         totalThick = ceil(max(max(stratCrossSec))) - floor(min(min(stratCrossSec))); % find total thickness of cross section in m
         cellThick = round(totalThick * res); % convert total thickness in cells
         timestepStrat = zeros(cellThick, domainDistance); % initiate timestep strat
         faciesStrat = zeros(cellThick,domainDistance); % initiate facies strat
         faciesCrossSec = totalFaciesSlice; 
         minElev=floor(min(stratCrossSec(:)));
         maxElev=ceil(max(stratCrossSec(:)));
         x=1:dx:(domainDistance*dx);
         y=(minElev+res_m/2):res_m:(maxElev-res_m/2);
         [ElevX, ElevY]=meshgrid(x,y);
     end

     % Initiate Starting Surface and Accumulator to fill Timestep and Facies Strat
     elevDiffCell = zeros(numSteps,domainDistance);
     lostToRound = zeros(1,domainDistance);
     accumulator = zeros(1,domainDistance);
     surfInit = stratCrossSec(stratstart,:); % starting surface and places in the correct pos %% [] should be 3? ** [] set to 4 %% DAE - totally confused why we choose 4?!
     numCols = length(surfInit); % total number of columns to loop through for next step 
     surfInitCells = fix(surfInit*res); % changes units from m to cells and rounds towards 0 
     lostToRoundSurf = surfInit - surfInitCells/res; % how many meters are lost to rounding in the second timestep
     currentSurf = surfInitCells; % initialize current surf as the initial surface in timestep 2
     currentSurf=(currentSurf-min(currentSurf))+10;
     accumulator = accumulator + lostToRoundSurf; % adds the amount lost to rounding in the second timestep to the accumulator 

    %% Apply acccumulator
     for i = stratstart:numSteps  % USUALLY STARTS AT 3!!!!!!!  []was at 2, I put it to 4
         elevDiffCell(i,:) = fix((elevDiff(i,:)*res)); % changes each row of total Elev Diff from m to cells (rounds towards 0)
         lostToRound = elevDiff(i,:) - elevDiffCell(i,:)/res; % the value of elevation in each cell lost to rounding during one timestep
         accumulator = accumulator + lostToRound; % apply lost to round to accumulator
         currentElevDiffCell = elevDiffCell(i,:); % change current row of elevDiffCell to a vector
         currentElevDiffCell(accumulator > 1/res) = currentElevDiffCell(accumulator > 1/res) + 1; % aggrade by one when accumulator is greater than 1 cell
         currentElevDiffCell(accumulator < (-1/res)) = currentElevDiffCell(accumulator < (-1/res)) - 1; % incise by one when accumulator is less than -1 cell
         elevDiffCell(i,:) = currentElevDiffCell; % applies changes back to elevDiffCell
         accumulator(accumulator < (-1/res)) = accumulator(accumulator < (-1/res)) + 1/res; % erase overfill
         accumulator(accumulator > 1/res) = accumulator(accumulator > 1/res) - 1/res; % erase overfill

         % Update strat with changes in cellular elevations
         for j = 1:domainDistance 
             if elevDiffCell(i,j) > 0 && currentSurf(j)+1 > 0 % if aggrading %% [] why is term 2 in this condition? **
                 timestepStrat(currentSurf(j)+1:(currentSurf(j)+elevDiffCell(i,j)),j) = i; % places values of timestep where aggradation is occuring
                 if faciesCrossSec(i,j) == 2 % if cell is active channel
                     faciesStrat(currentSurf(j)+1:currentSurf(j) + elevDiffCell(i,j),j) = 3;

                 elseif faciesCrossSec(i,j) == 1 % if cell is abandoned
                     faciesStrat(currentSurf(j)+1:currentSurf(j) + elevDiffCell(i,j),j) = 2;

                 elseif faciesCrossSec(i,j) == 0 % if cell is "floodplain and healed"
                     faciesStrat(currentSurf(j)+1:currentSurf(j) + elevDiffCell(i,j),j) = 1;
                 end
             elseif elevDiffCell(i,j) < 0 && currentSurf(j) + elevDiffCell(i,j)+1 > 0 % if incising %% [] why is term 2 in this condition? **
                 timestepStrat((currentSurf(j) + elevDiffCell(i,j)+1):currentSurf(j),j) = 0; % erases correlating cells with 0 or where incision is occuring %%DAE - WHY the "+1"
                 faciesStrat(currentSurf(j) + elevDiffCell(i,j)+1:currentSurf(j),j) = 0; % should do the same to the facies
             else
             end
              currentSurf(:,j) = currentSurf(:,j) + elevDiffCell(i,j); % updates the current surface after each j iteration
         end
     end

    %% Save Facies, Timestep, and totalelevdiff Matrices
    % Need to figure out how to not hardcode in the file name later
    myfilename_1 = strcat('faciesStrat',num2str(sliceLocal),'.mat');
    fullfilename_1 = fullfile(exportDirectory,myfilename_1);
    save(fullfilename_1,'faciesStrat', 'x', 'y');

    myfilename_2 = strcat('timestepStrat',num2str(sliceLocal),'.mat');
    fullfilename_2 = fullfile(exportDirectory,myfilename_2);
    save(fullfilename_2,'timestepStrat');

    myfilename_3 = strcat('totalFacies',num2str(sliceLocal),'.mat');
    fullfilename_3 = fullfile(exportDirectory,myfilename_3);
    save(fullfilename_3,'totalFaciesSlice');


    %% Modal Percentages for Facies
%     if runNum == 1
%         bedloadf_perc = zeros(1,totalRunNum); % intitialize array to collect bl values for cross sections
%         overbankf_perc = zeros(1,totalRunNum); % intitialize array to collect bl values for cross sections
%         floodplain_perc = zeros(1,totalRunNum); % intitialize array to collect bl values for cross sections
%         occ_blf = zeros(1,totalRunNum);
%         occ_obf = zeros(1,totalRunNum);
%         occ_fp = zeros(1,totalRunNum);
%         tot_fac = zeros(1,totalRunNum);
%     end
%     occ_blf(runNum) = length(find(faciesStrat == 3)); % how many occurences of bedload fill facies there are in the cross section
%     occ_obf(runNum) = length(find(faciesStrat == 2)); % how many occurences of overbank fill facies there are the cross section
%     occ_fp(runNum) = length(find(faciesStrat == 1)); % how many occurences of floodplain facies there are the cross section
%     tot_fac(runNum) = occ_blf(runNum) + occ_obf(runNum) + occ_fp(runNum); % total the amount of facies 
%     bedloadf_perc(runNum) = (occ_blf(runNum)/tot_fac(runNum))*100;
%     overbankf_perc(runNum) = (occ_blf(runNum)/tot_fac(runNum))*100;
%     floodplain_perc(runNum) = (occ_blf(runNum)/tot_fac(runNum))*100;
%     bedloadf_perc(runNum) = (occ_blf(runNum)/tot_fac(runNum))*100;
%     overbankf_perc(runNum) = (occ_obf(runNum)/tot_fac(runNum))*100;
%     floodplain_perc(runNum) = (occ_fp(runNum)/tot_fac(runNum))*100;
    
    %% Model Output Plots

% %     if toggle == 0     
% %         figure(1)
% %         fig1 = imagesc(timestepStrat); % x is width and y is cells
% %         set(gca,'YDir','normal')
% %         xlabel('Distance (units of 500m)')
% %         ylabel('#Cells')
% %         title(strcat('Timestep Stratigraphy X-Section:', num2str(sliceLocal)));
% %         a = colorbar;
% %         a.Label.String = 'Timestep';
% % %         saveas(fig1,fullfile(exportDirectory2,['timestepStrat_new' num2str(exportName) '_Res' num2str(res) '_' datestr(now,'HH-MM-SS') '.fig']));
% % %         saveas(fig1,fullfile(exportDirectory2,['timestepStrat_new' num2str(exportName) '_Res' num2str(res) '_' datestr(now,'HH-MM-SS') '.png']));
% % 
% %         % Plot faciesStrat
        figure(3)
        fig3 = imagesc(x,y,faciesStrat);
        hold on 
        set(gca,'YDir','normal')
        xlabel('Distance (units of 500m)')
        ylabel('Thickness (m)')
        title(strcat('Facies Stratigraphy X-Section:', num2str(sliceLocal)));
        colormap(faciesColors)
        a = colorbar;
        a.Label.String = 'Facies';
% % %         saveas(fig3,fullfile(exportDirectory2,['faciesStrat_new' num2str(exportName) '_Res' num2str(res) '_' datestr(now,'HH-MM-SS') '.fig']));
% % %         saveas(fig3,fullfile(exportDirectory2,['faciesStrat_new' num2str(exportName) '_Res' num2str(res) '_' datestr(now,'HH-MM-SS') '.png']));
        
        %%end
    end
% %% Plot modal percentages
% figure(7)
%     fig7 = plot(sliceLocal_input(1:runNum),bedloadf_perc, 'g.-');
%     hold on
%     plot(sliceLocal_input(1:runNum), overbankf_perc, 'b.-')
%     plot(sliceLocal_input(1:runNum),floodplain_perc,'r.-')
%     hold off 
%     set(gca,'Xdir','reverse')
%     title('Model Run Modal Percentages')
%     xlabel('Cross Section Location')
%     ylabel('Modal Percentage')
%     legend('Bedload Fill', 'Overbank Fill','Floodplain')
%     grid on
%     saveas(fig7,fullfile(exportDirectory,['modalperc' num2str(exportName) '.fig']));
%     saveas(fig7,fullfile(exportDirectory,['modalperc' num2str(exportName) '.epsc']));
%     
% figure(8)
%     fig8 = plot(sliceLocal_input(1:runNum),overbankf_perc,'b.-');
%     set(gca,'Xdir','reverse')
%     title('Overbank Modal Percentage')
%     xlabel('Cross Section Location')
%     ylabel('Modal Percentage')
%     grid on 
%     saveas(fig8,fullfile(exportDirectory,['modalperc_overbank' num2str(exportName) '.fig']));
%     saveas(fig8,fullfile(exportDirectory,['modalperc_overbank' num2str(exportName) '.epsc']));
%     
%     close all
    
% %% Save Modal Percentages for Run   
% myfilename_4 = strcat('BedloadPerc.mat');
% fullfilename_4 = fullfile(exportDirectory,myfilename_4);
% save(fullfilename_4,'bedloadf_perc');
% 
% myfilename_5 = strcat('OverbankPerc.mat');
% fullfilename_5 = fullfile(exportDirectory,myfilename_5);
% save(fullfilename_5,'overbankf_perc');
% 
% myfilename_6 = strcat('FloodplainPerc.mat');
% fullfilename_6 = fullfile(exportDirectory,myfilename_6);
% save(fullfilename_6,'floodplain_perc');

end
%% Timer
% % toc
% % 
% %         Topo = thisStrike'; % flip to input into function
% %         blockSize = 100; % add slices together at every 100 timesteps 
% %         numSteps = length(Topo(1,:));
% %         Topo_new = zeros(length(Topo(:,1)),floor((numSteps/blockSize))); % initialize to resize together at every 100 timesteps
% %         n = length(Topo_new(1,:));
% % %         clear newSurf
% % 
% %         %% Resizing Loop
% %         % loop to resize input topographic elevations through time to stitch together ever 100th layer thru the whole volume of the fan: makes computation for compensation metric much more efficient
% %         for i = 1:n-1
% %             m = (1:blockSize:numSteps);
% %             Topo_new(:,i) = Topo(:,m(i)); %+ Topo(:,m(i+1));
% %         end
% %         clear Topo
% %         Topo_new = Topo_new(:,1:end-1); % remove last empty surface
% %         
% %         % fix newSurf so elevation values only increase through time, takes about 2.5 hr to run this loop with about  
% %         % sort columns for every row   
% %         for r = 1:length(Topo_new(:,1)) 
% %             for c = 1:length(Topo_new(1,:))
% %                 Topo_new(r,c) = min(Topo_new(r,c:end)); % for each column this loop places the lowest elevation values in the row at the lowest column, essentially sorts elevation data
% %             end
% %         end
% %         
% %         hold on
% %         %figure(4)
% %         plot(x,Topo_new(:,1:10:end),'k'); % visualize topographic surfaces
% %         ylabel('Thickness (m)')
% %         xlabel('Distance (units of 500m)')
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        