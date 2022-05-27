function animalSuffix = getAnimalSuffix(animalCodes)
% only apply to 5CSRTT
% Get animal suffix eg. '_123A'
allAnimals = {'0171','0179','0180','0181'};
mask = ismember(allAnimals,animalCodes);
temp = '1234';
id = temp(mask);
animalSuffix = ['_' id 'A'];
end