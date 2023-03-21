function LH = likelihood(choice, probChoice)
% computes likelihood for learning models
if size(choice,1)~=size(probChoice,1)
    if length(choice) == length(probChoice)
        choice = choice';
    end
end
chosenProb = (probChoice.^(choice)).*((1-probChoice).^(1-choice));
chosenProb(chosenProb==0)= 0.0001;
LH = -sum(log(chosenProb));

end