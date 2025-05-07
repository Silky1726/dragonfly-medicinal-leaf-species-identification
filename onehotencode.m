% Assuming 'score' is a matrix of probabilities or raw scores
[numClasses, numSamples] = size(score);
oneHotEncodedResult = zeros(numClasses, numSamples);

for i = 1:numSamples
    [~, classLabel] = max(score(:, i));
    oneHotEncodedResult(classLabel, i) = classLabel;
end
