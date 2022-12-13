function state = observatorFnc(observation_data, control_data, time_samples, G1_vec, G2_vec)
    % funkcja licząca numerycznie obserwacje
    % observation_data - zebrane obserwacje
    % control_data - zebrane_sterowania
    % time_samples - próbki czasu, dla których zostały wyliczone macierze G1 G2
    % G1_vec - macierz G1
    % G2_vec - macierz G2
    
    G1xY = squeeze(bsxfun(@times,G1_vec,reshape(observation_data, 1, 1, [])));
    G2xU = squeeze(bsxfun(@times,G2_vec,reshape(control_data, 1, 1, [])));
    
    state_obs = trapz(time_samples, G1xY, 2);
    control_obs = trapz(time_samples, G2xU, 2);
    
    state = state_obs + control_obs;
end
