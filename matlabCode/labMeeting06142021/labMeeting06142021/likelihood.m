function LH = likelihood(choice, probChoice)
% computes likelihood for learning models

chosenProb = (probChoice.^(choice)).*((1-probChoice).^(1-choice));
chosenProb(chosenProb==0)= 0.0001;
LH = -sum(log(chosenProb));

end