function fitStanModels(varargin)

p = inputParser;
% default parameters if none given
p.addParameter('animals', 'ninety')
p.addParameter('modelName', 'sixParam_absPePeAN_bi')
p.addParameter('params', {'aNscale', 'aNmin', 'aP', 'aF', 'aPE', 'beta'})
p.addParameter('lickBeh', 'c');
p.addParameter('bernFlag', 1)
p.addParameter('tWflag', 0)
p.addParameter('data', [])
p.addParameter('iter', 10000)
p.addParameter('beh', 'clean')
p.parse(varargin{:});

pPflag = 0;

if iscell(p.Results.animals)
    animals = p.Results.animals;
    pPflag = 1;
else
    switch p.Results.animals
        case 'ninety'
            animals = [{'CG45' 'CG46' 'CG47' 'CG48' 'CG49' 'CG50' 'CG51' 'CG52' 'CG53' 'CG54' ...
                        'CG55' 'CG56' 'CG57' 'CG58' 'CG59' 'CG60' 'CG61' 'CG62' 'CG63' 'CG64' ...
                        'CG65' 'CG66' 'CG67' 'CG68' 'CG70' 'CG75' 'CG77' 'CG78' 'CG79' 'CG80' ... 
                        'CG81' 'CG82' 'CG83' 'CG84' 'CG85' 'CG86' 'CG87' 'CG88' 'CG89' 'CG90' ... 
                        'CG91' 'CG92' 'CG93' 'CG94' 'CG95' 'TP02' 'TP03'}];
        case 'meow'
            animals = [{'CG78' 'CG79' 'CG80' ... 
                        'CG81' 'CG82' 'CG83' 'CG84' 'CG85' 'CG86' 'CG87' 'CG88' 'CG89' 'CG90' ... 
                        'CG91' 'CG92' 'CG93' 'CG94' 'CG95' 'TP02' 'TP03'}];
        case 'ninetyNoRec'
            animals = [{'CG45' 'CG46' 'CG47' 'CG48' 'CG49' 'CG50' 'CG52' 'CG53' 'CG54' 'CG55' ...
                        'CG56' 'CG57' 'CG58' 'CG67' 'CG68' 'CG70' 'CG77' 'CG78' 'CG80' 'CG82' ... 
                        'CG83' 'CG84' 'CG85' 'CG86' 'CG87' 'CG88' 'CG89'}];
        case 'recordings'
            animals = [{'CG09' 'CG16' 'CG75' 'CG79' 'CG81'}];
        case 'opto'
            animals = [{'CG93' 'CG94' 'CG95' 'TP02' 'TP03'}];
        case 'pP_5ht'
            animals = [{'BB025', 'BB031', 'BB033', 'BB035', 'BB036'}];
            pPflag = 1;
        case 'pP_3probs'
            animals = [{'SL01' ,'SL02', 'SL03', 'SL04'}];
            pPflag = 1;
    end
end


if pPflag == 1 & isempty(p.Results.data)
    load('C:\Users\cooper\Desktop\pP\fullStruct.mat')
elseif pPflag == 1 & ~isempty(p.Results.data)
    fullStruct = p.Results.data;
end    
    
mdlTxt = [];
if pPflag
    for i = 1:length(animals)
      tmpStr = ['stan_qLearningFit_pP(''%s'', ''iter'', %d, ''modelName'', ''%s'', ''data'', fullStruct, ''lickBeh'', ''%s'', ''paramNames'', {' ...
          [repmat(' ''%s'', ', 1, length(p.Results.params))] '});'];
      tmp = sprintf(tmpStr, animals{i}, p.Results.iter, p.Results.modelName, p.Results.lickBeh, p.Results.params{:});
      mdlTxt = [mdlTxt tmp];
    end
elseif p.Results.tWflag
    for i = 1:length(animals)
      tmpStr = ['stan_qLearningFit_tranWeight(''goodBehDays.xlsx'', ''%s'', ''%s'', ''iter'', %d, ''bernFlag'', %d, ''modelName'', ''%s'', ''paramNames'', {' ...
          [repmat(' ''%s'', ', 1, length(p.Results.params))] '});'];
      tmp = sprintf(tmpStr, animals{i}, p.Results.beh, p.Results.iter, p.Results.bernFlag, p.Results.modelName, p.Results.params{:});
      mdlTxt = [mdlTxt tmp];
    end
else
    for i = 1:length(animals)
      tmpStr = ['stan_qLearningFit(''goodBehDays.xlsx'', ''%s'', ''%s'', ''iter'', %d, ''bernFlag'', %d, ''modelName'', ''%s'', ''paramNames'', {' ...
          [repmat(' ''%s'', ', 1, length(p.Results.params))] '});'];
      tmp = sprintf(tmpStr, animals{i}, p.Results.beh, p.Results.iter, p.Results.bernFlag, p.Results.modelName, p.Results.params{:});
      mdlTxt = [mdlTxt tmp];
    end
end

eval(mdlTxt);

end

