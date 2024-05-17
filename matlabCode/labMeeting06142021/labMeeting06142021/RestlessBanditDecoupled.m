classdef RestlessBanditDecoupled < handle
    %
    %   RestlessBanditDecoupled class
    %
    % A class that simulates the VR-VR dynamic foraging task
    % Initialization: Takes as inputs reward probabilities, block lengths, and number of trials
    % Running: Takes as inputs the choice
    %
    % BAB 6/21/17
    % modified by CG 5/29/18
    %
    properties
        RewardProbabilitiesList
        RewardProbabilities % left, right
        BlockLength
        BlockEndL
        BlockEndR
        AllRewards = [];
        AllChoices = [];
        BlockSwitch_Flag = true;
        NewBlockL_Flag = false;
        NewBlockR_Flag = false;
        BlockSwitchL = [];
        BlockSwitchR = [];
        BlockProbs = [];
        Trial
        MaxTrials
        RandomSeed
        BlockStagger
        HigherCount
        PersevCount
        ITI
    end
    methods
        function obj = RestlessBanditDecoupled(varargin) % constructor            
            p = inputParser;
            
            % default parameters if none given
            p.addParameter('RewardProbabilities', [90 50 10]);
            p.addParameter('BlockLength', [20 35]);
            p.addParameter('MaxTrials', 1000);
            p.addParameter('RandomSeed', 1);
            p.parse(varargin{:});
            
            obj.RandomSeed = p.Results.RandomSeed;
            rng(obj.RandomSeed);
            obj.RewardProbabilitiesList = p.Results.RewardProbabilities;
            obj.RewardProbabilities(1,1) = obj.RewardProbabilitiesList(randi(3));
            obj.RewardProbabilities(1,2) = obj.RewardProbabilitiesList(randi(3));
            while obj.RewardProbabilities(1,:) == 10
                obj.RewardProbabilities(1,1) = obj.RewardProbabilitiesList(randi(3));
                obj.RewardProbabilities(1,2) = obj.RewardProbabilitiesList(randi(3));
            end
            obj.BlockProbs = obj.RewardProbabilities;
            obj.BlockLength = p.Results.BlockLength;
            obj.MaxTrials = p.Results.MaxTrials;
            obj.AllRewards = NaN(obj.MaxTrials, 2);
            obj.AllChoices = NaN(obj.MaxTrials, 2); % indexed as left and right
            obj.Trial = 0;
            obj.BlockSwitchL = 1;
            obj.BlockSwitchR = 1;
            if obj.RewardProbabilities(1,1) > obj.RewardProbabilities(1,2)
                obj.HigherCount = [1 0];
            elseif obj.RewardProbabilities(1,1) < obj.RewardProbabilities(1,2)
                obj.HigherCount = [0 1];
            else
                obj.HigherCount = [1 1];
            end
            obj.BlockEndL = randi([min(obj.BlockLength) max(obj.BlockLength)]) + obj.Trial;
            obj.BlockEndR = randi([min(obj.BlockLength) max(obj.BlockLength)]) + obj.Trial;
            obj.BlockSwitchL = [obj.BlockSwitchL obj.BlockEndL];
            obj.BlockSwitchR = [obj.BlockSwitchR obj.BlockEndR];
            obj.BlockStagger = ceil(mean(p.Results.BlockLength) / 2);
            if binornd(1,0.5) == 0
                obj.BlockEndL = obj.BlockEndL - obj.BlockStagger;
                obj.BlockSwitchL(end) = obj.BlockEndL;
            else
                obj.BlockEndR = obj.BlockEndR - obj.BlockStagger;
                obj.BlockSwitchR(end) = obj.BlockEndR;
            end
            obj.PersevCount = [0 0];
            obj.ITI = 0;
        end
        function obj = inputChoice(obj, currChoice)
            % increment trial
            obj.Trial = obj.Trial + 1;
            
            % turn off flag used for plotting block switches
            if obj.BlockSwitch_Flag == true & obj.Trial > 1
                obj.BlockSwitch_Flag = false;
            end 
            
            % persevAdd code
            if currChoice == [1 0] & obj.RewardProbabilities(1,1) == 10
                obj.PersevCount(1,:) = [obj.PersevCount(1,1) + 1 0];
            elseif currChoice == [0 1] & obj.RewardProbabilities(1,2) == 10
                obj.PersevCount(1,:) = [0 obj.PersevCount(1,2) + 1];
            end
            if obj.PersevCount(1,1) >= 4 & obj.RewardProbabilities(1,1) == 10
                obj.BlockEndL = obj.BlockEndL + obj.PersevCount(1,1);
                obj.BlockSwitchL(end) = obj.BlockEndL;
                obj.BlockEndR = obj.BlockEndR + obj.PersevCount(1,1);
                obj.BlockSwitchR(end) = obj.BlockEndR;
                obj.PersevCount(1,1) = 0;
                obj.NewBlockL_Flag = false;
                obj.NewBlockR_Flag = false;
            elseif obj.PersevCount(1,2) >= 4 & obj.RewardProbabilities(1,2) == 10
                obj.BlockEndL = obj.BlockEndL + obj.PersevCount(1,2);
                obj.BlockSwitchL(end) = obj.BlockEndL;
                obj.BlockEndR = obj.BlockEndR + obj.PersevCount(1,2);
                obj.BlockSwitchR(end) = obj.BlockEndR;
                obj.PersevCount(1,2) = 0;
                obj.NewBlockL_Flag = false;
                obj.NewBlockR_Flag = false;
            end
                
            % currChoice needs to be a 1x2 binary vector for [L R] choice
            if obj.NewBlockL_Flag == true
                obj = generateBlockL(obj);
            end
            
            if obj.NewBlockR_Flag == true
                obj = generateBlockR(obj);
            end
            
            % input choice
            Rwd_rand = randi([0 99]);
            
            if all(currChoice == [1 0]) % left choice
                obj.AllChoices(obj.Trial, :) = [1 0];
                if obj.RewardProbabilities(1) > Rwd_rand % harvest L reward
                    obj.AllRewards(obj.Trial, :) = [1 0];
                else
                    obj.AllRewards(obj.Trial, :) = [0 0];
                end
            elseif all(currChoice == [0 1]) % right choice
                obj.AllChoices(obj.Trial, :) = [0 1];
                if obj.RewardProbabilities(2) > Rwd_rand % harvest R reward
                    obj.AllRewards(obj.Trial, :) = [0 1];
                else
                    obj.AllRewards(obj.Trial, :) = [0 0];
                end
            else
                error('Choice should be [1 0] (left choice) or [0 1] (right choice)')
            end
            
            % next block
            if obj.Trial == obj.BlockSwitchL(end)
                obj.NewBlockL_Flag = true;
            end
            if obj.Trial == obj.BlockSwitchR(end)
                obj.NewBlockR_Flag = true;
            end
            
            obj.ITI = -log(rand(1))/0.3;
            
        end
        
        function obj = generateBlockL(obj)
            obj.NewBlockL_Flag = false;
            obj.BlockSwitch_Flag = true;
            obj.BlockEndL = randi([min(obj.BlockLength) max(obj.BlockLength)]) + obj.Trial;
            obj.BlockSwitchL = [obj.BlockSwitchL obj.BlockEndL];
            if obj.HigherCount(1) == 3
                obj.RewardProbabilities(1,1) = 10;
                obj.HigherCount(1) = 0;
                obj.HigherCount(2) = obj.HigherCount(2) + 1;
            else
                tmp = obj.RewardProbabilities(1,1);
                while obj.RewardProbabilities(1,1) == tmp;
                    obj.RewardProbabilities(1,1) = obj.RewardProbabilitiesList(randi(3));
                end
                if obj.RewardProbabilities(1,1) > obj.RewardProbabilities(1,2)
                    obj.HigherCount(1) = obj.HigherCount(1) + 1;
                    obj.HigherCount(2) = 0;
                elseif obj.RewardProbabilities(1,1) < obj.RewardProbabilities(1,2)
                    obj.HigherCount(2) = obj.HigherCount(2) + 1;
                    obj.HigherCount(1) = 0;
                else
                    obj.HigherCount(1) = obj.HigherCount(1) + 1;
                    obj.HigherCount(2) = obj.HigherCount(2) + 1;
                end
            end
            if obj.RewardProbabilities(1,:) == 10
                obj.BlockEndL = obj.BlockEndL - obj.BlockStagger;
                obj.BlockSwitchL(end) = obj.BlockEndL;
                obj.HigherCount(:) = 0;
                obj.BlockEndR = obj.Trial - 1;
                obj.BlockSwitchR(end) = obj.BlockEndR;
                obj = generateBlockR(obj);
            elseif obj.NewBlockR_Flag == false
                obj.BlockProbs = [obj.BlockProbs; obj.RewardProbabilities];
            end
        end

        function obj = generateBlockR(obj)
            obj.NewBlockR_Flag = false;
            obj.BlockSwitch_Flag = true;
            obj.BlockEndR = randi([min(obj.BlockLength) max(obj.BlockLength)]) + obj.Trial;
            obj.BlockSwitchR = [obj.BlockSwitchR obj.BlockEndR];
            if obj.HigherCount(2) == 3
                obj.RewardProbabilities(1,2) = 10;
                obj.HigherCount(2) = 0;
                obj.HigherCount(1) = obj.HigherCount(1) + 1;
            else
                tmp = obj.RewardProbabilities(1,2);
                while obj.RewardProbabilities(1,2) == tmp;
                    obj.RewardProbabilities(1,2) = obj.RewardProbabilitiesList(randi(3));
                end
                if obj.RewardProbabilities(1,2) > obj.RewardProbabilities(1,1)
                    obj.HigherCount(2) = obj.HigherCount(2) + 1;
                    obj.HigherCount(1) = 0;
                elseif obj.RewardProbabilities(1,2) < obj.RewardProbabilities(1,1)
                    obj.HigherCount(1) = obj.HigherCount(1) + 1;
                    obj.HigherCount(2) = 0;
                else
                    obj.HigherCount(2) = obj.HigherCount(2) + 1;
                    obj.HigherCount(1) = obj.HigherCount(1) + 1;
                end
            end
            if obj.RewardProbabilities(1,:) == 10
                obj.BlockEndR = obj.BlockEndR - obj.BlockStagger;
                obj.BlockSwitchR(end) = obj.BlockEndR;
                obj.HigherCount(:) = 0;
                obj.BlockEndL = obj.Trial - 1;
                obj.BlockSwitchL(end) = obj.BlockEndL;
                obj = generateBlockL(obj);
            else
                obj.BlockProbs = [obj.BlockProbs; obj.RewardProbabilities];
            end
        end
    end
end



        