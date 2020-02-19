function LH = likelihood(choice, probChoice)
% computes likelihood for learning models


Pf = choice.* probChoice;
nonzeros = Pf ~= 0;

if sum(sum(nonzeros)) < size(choice, 1) % if there are probabilities = 0
    LH = -1 * sum(log(0.01*ones(10*size(choice, 1), 1))); % give each choice a probability of 1% to give a bad likelihood
else
    Pf = Pf(nonzeros);
    LH = -1 * sum(log(Pf));
end