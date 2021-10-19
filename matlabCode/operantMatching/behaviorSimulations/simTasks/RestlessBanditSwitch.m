classdef RestlessBanditSwitch < handle
    %
    %   RestlessBandit class
    %
    % A class that simulates the VR-VR dynamic foraging task
    % Initialization: Takes as inputs reward probabilities, block lengths, and number of trials
    % Running: Takes as inputs the choice
    %
    % BAB 6/21/17
    % modified by CG 3/22/18
    %
    properties
        RewardProbabilities % left, right
        RewardProbabilitiesListPre = [];
        RewardProbabilitiesListPost = [];
        BlockLength
        BlockEnd
        AllRewards = [];
        AllChoices = [];
        BlockSwitch_Flag = true;
        NewBlock_Flag = false;
        BlockSwitch = [];
        BlockInd
        BlockCount
        BlockProbs = [];
        Trial
        MaxTrials
        RandomSeed
        PersevCount
    end
    methods
        function obj = RestlessBanditSwitch(varargin) % constructor            
            p = inputParser;
            
            % default parameters if none given
            p.addParameter('RewardProbabilitiesListPre', [90 10; 10 90]);
            p.addParameter('RewardProbabilitiesListPost', [50 10; 10 50]);
            p.addParameter('BlockLength', [20 35]);
            p.addParameter('MaxTrials', 1000);
            p.addParameter('RandomSeed', 1);
            p.parse(varargin{:});
            
            obj.RandomSeed = p.Results.RandomSeed;
            rng(obj.RandomSeed);
            obj.RewardProbabilitiesListPre = p.Results.RewardProbabilitiesListPre;
            obj.RewardProbabilities = p.Results.RewardProbabilitiesListPre(1,:);
            obj.RewardProbabilitiesListPost = p.Results.RewardProbabilitiesListPost;
            obj.BlockProbs = obj.RewardProbabilities;
            obj.BlockLength = p.Results.BlockLength;
            obj.MaxTrials = p.Results.MaxTrials;
            obj.AllRewards = NaN(obj.MaxTrials, 2);
            obj.AllChoices = NaN(obj.MaxTrials, 2); % indexed as left and right
            obj.Trial = 0;
            obj.BlockSwitch = 1;            
            obj.BlockEnd = randi([min(obj.BlockLength) max(obj.BlockLength)]) + obj.Trial;
            obj.BlockSwitch = [obj.BlockSwitch obj.BlockEnd];
            obj.BlockInd = 1;
            obj.BlockCount = 1;
            obj.PersevCount = [0 0];
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
                obj.BlockEnd = obj.BlockEnd + obj.PersevCount(1,1);
                obj.BlockSwitch(end) = obj.BlockEnd;
                obj.PersevCount(1,1) = 0;
                obj.NewBlock_Flag = false;
            elseif obj.PersevCount(1,2) >= 4 & obj.RewardProbabilities(1,2) == 10
                obj.BlockEnd = obj.BlockEnd + obj.PersevCount(1,1);
                obj.BlockSwitch(end) = obj.BlockEnd;
                obj.PersevCount(1,2) = 0;
                obj.NewBlock_Flag = false;
            end
            
            % generate a new block
            if obj.NewBlock_Flag == true
                obj.BlockEnd = randi([min(obj.BlockLength) max(obj.BlockLength)]) + obj.Trial;
                obj.BlockSwitch = [obj.BlockSwitch obj.BlockEnd];
                obj.BlockInd = obj.BlockInd + 1;
                if obj.BlockInd > length(obj.RewardProbabilitiesListPre)
                    obj.BlockInd = 1;
                end
                obj.BlockCount = obj.BlockCount + 1;
                if obj.BlockCount < 3
                    obj.RewardProbabilities = obj.RewardProbabilitiesListPre(obj.BlockInd,:); % switch reward probabilities
                else
                    obj.RewardProbabilities = obj.RewardProbabilitiesListPost(obj.BlockInd,:);
                end
                obj.BlockProbs = [obj.BlockProbs; obj.RewardProbabilities];
                obj.NewBlock_Flag = false;
                obj.BlockSwitch_Flag = true;
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
            if obj.Trial == obj.BlockSwitch(end)
                obj.NewBlock_Flag = true;
            end
        end
    end
end