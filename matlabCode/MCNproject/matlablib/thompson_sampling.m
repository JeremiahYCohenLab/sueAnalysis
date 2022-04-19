function choice = thompson_sampling(a1, b1, a2, b2)
    beta1 = betarnd(a1, b1);
    beta2 = betarnd(a2, b2);

    choice = double(beta2 > beta1);