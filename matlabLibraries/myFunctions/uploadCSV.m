function uploadCSV(sessionList)

dataDir = '/allen/aind/scratch/ephys/temp/Behavior/';
localDir = 'Z:\ephys\temp\Behavior';

platform = cell(length(sessionList), 1);
acq_datetime = cell(length(sessionList), 1);
subject_id = cell(length(sessionList), 1);
s3_bucket = cell(length(sessionList), 1);
modality0 = cell(length(sessionList), 1);
source0 = cell(length(sessionList), 1);
modality1 = cell(length(sessionList), 1);
source1 = cell(length(sessionList), 1);
modality2 = cell(length(sessionList), 1);
source2 = cell(length(sessionList), 1);


for i = 1:length(sessionList)
    session = sessionList{i};
    animalName = session(1:6);
    date = session(8:17);
    time = session(19:26);
    
    % conver time format
    time24 = [time(1:2) ':' time(4:5) ':' time(7:8)];

    platform{i} = 'behavior';
    acq_datetime{i} = [date ' ' time24];
    subject_id{i} = animalName;
    s3_bucket{i} = 'aind-behavior-data';

    modality0{i} = 'behavior-videos';
    source0{i} = [dataDir animalName '/' session '/VideoFolder'];

    modality1{i} = 'trained-behavior';
    source1{i} = [dataDir animalName '/' session '/BehaviorFolder'];

    modality2{i} = 'ecephys';
    source2{i} = [dataDir animalName '/' session '/EphysFolder/'];
    allFiles = dir([localDir '\' animalName '\' session '\EphysFolder\']);
    ephysFolder = {allFiles.name};
    ephysFolder = ephysFolder{contains(ephysFolder, animalName)};
    source2{i} = [source2{i} ephysFolder];


    % merge behavior folder

    if ~exist([localDir '\' animalName '\' session '\BehaviorFolder'], 'dir')
        mkdir([localDir '\' animalName '\' session '\BehaviorFolder']);
        movefile([localDir '\' animalName '\' session '\TrainingFolder'], [localDir '\' animalName '\' session '\BehaviorFolder']);
        movefile([localDir '\' animalName '\' session '\HarpFolder'], [localDir '\' animalName '\' session '\BehaviorFolder']);
    end

end

tbl = table(platform, acq_datetime, subject_id, s3_bucket, modality0, source0, modality1, source1, modality2, source2);

tbl = renamevars(tbl,{'source0', 'source1', 'source2'}, {'modality0.source', 'modality1.source', 'modality2.source'});

writetable(tbl, 'F:\uploadData.csv');
