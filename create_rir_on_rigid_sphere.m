% This script creates room impulse response from sound source to mics on a rigid ball in self-defined conditions.
% num: number of RIRs you want to generate
% ===================================================
% !NOTE: Please check all the parameters and change them properly before you use them
% ===================================================
clc
clear
close

PARA = ~false;
CREATE_LIB = false;
DEBUG = ~true;
if DEBUG
    fprintf('!Debug mode is activated');
end
if(PARA)
    cc = parcluster('local');
    ncores = feature('numCores');
    cc.NumWorkers = ncores;
    parpool(cc, cc.NumWorkers);
end

% Number of TFs to be generated
tf_num = 936;
% Note: for sample balance, each position has the same number of TFs
% For 72 azimuths and 13 elevations, there are 72 * 13 = 936 positions in
% total. If we generate 500 samples (rooms) for each position, then there
% would be 936 * 500 = 468000 Tfs in total.
% In the following implementaion, the certain position will be determined
% by the index of parfor loop.

dist_in_lib = [0.5:0.05:0.95, 1:0.1:3, 3:0.5:10, 10:1:20]'; % These distances would be considered when building the TF lib
reso_in_deg = 5;
resolution = reso_in_deg / 180 * pi;
angle_in_lib = [0:resolution:pi]';
fs = 8e3;         
tf_len = 4096; 
c = 344; % m/s
nMic = 32;
nAzi = 72;
nEle = 13;
total_pos = nAzi * nEle;
% reflection order
order = 3;

aziList = linspace(0, 355, nAzi) / 180 * pi;
eleList = linspace(-30, 30, nEle) / 180 * pi;

r = 0.042; % The radius of eigenMike microphone array
array_z_range   = [0.8, 1.5]; %(m), array height
array_xy_jittering = 0.6; %(m)
T60_range       = [0.2, 0.6]; %s

room_l_range    = [6, 10]; %(m)
room_w_range    = [6, 10]; %(m)
room_h_range    = [2.5, 4];  %(m)
% Thresholds
array_dist_thres      = 0.5;
source_dist_thres = 0.5;
source_z_thres = 0.3;
source_center_range = [1, 3];

MicCoor_raw = mic_angle(); % [azimuth, elevation],azimuth[0~2pi], elevation[-pi/2, pi/2]
% MicCoor = load('mic.mat');
MicCoorRel = [MicCoor_raw, repelem(r, nMic)'];   % spherical coordinates of mics relative to sphere center
MicCoorRel = sph2cart(MicCoorRel(:, 1), MicCoorRel(:, 2), MicCoorRel(:, 3));


if(CREATE_LIB)
    TF_Lib = CalculateTF_Lib(fs, tf_len, angle_in_lib, dist_in_lib);
    % calculate the TF LIBRARY under free sound field
    save(['eigenMike_TF_Lib_4096_fs_', num2str(fs), '.mat'],  'TF_Lib');
end

TF_Lib = load(['eigenMike_TF_Lib_4096_fs_', num2str(fs), '.mat']);
TF_Lib = TF_Lib.TF_Lib;
rir_path = '/mnt/Disk3/cjf/SpeechEnhancement/EarlyReflections/210207eigenMikeTF/tr/';
if(~isfolder(rir_path));mkdir(rir_path);end

% !Note: If you are in debug mode, the following loop should be for,
% instead of parfor.

parfor room_idx = 0:tf_num-1
    room_idx
    T60             = draw(T60_range);  % Reverberation time (s)
    % Generate the locations of rooms, mic array and sound source randomly until conditions are satisfied
   
    pos_base = floor(room_idx / total_pos);
    pos_idx = mod(room_idx, total_pos);
    row = floor(pos_idx / nEle);
    col = mod(pos_idx, nEle);
    
    azi = aziList(row + 1);
    ele = eleList(col + 1);
    while 1
        room_l  = draw(room_l_range);   
        room_w  = draw(room_w_range);
        room_h  = draw(room_h_range);

        array_x = room_l/2 + (rand*2-1)*array_xy_jittering; % array_x is +-array_xy_jittering near the center
        array_y = room_w/2 + (rand*2-1)*array_xy_jittering; % array_y is +-array_xy_jittering near the center
        array_z = draw(array_z_range);
        
        if array_x <= array_dist_                                                                     thres || array_x >= room_l - array_dist_thres || array_y <= array_dist_thres || array_y >= room_w - array_dist_thres
            continue
        end

        source_center_dist = draw(source_center_range);
       
        source_image = source_center_dist * cos(ele);
        source_position_x = source_image * cos(azi) + array_x;
        source_position_y = source_image * sin(azi) + array_y;
        source_position_z = source_center_dist * sin(ele) + array_z;   %

        if source_position_x >= source_dist_thres && source_position_x <= room_l - source_dist_thres && source_position_y >= source_dist_thres && source_position_y <= room_w - source_dist_thres && source_position_z >= source_z_thres && source_position_z <= room_h - source_z_thres
            break
        end
        if DEBUG
            fprintf('room x:%.3fm, y:%.3fm, z:%.3fm\n', room_l, room_w, room_h);
            fprintf('source x:%.3fm, y:%.3fm, z:%.3fm\n', source_position_x, source_position_y, source_position_z);
        end
    end
    disp('Condition chosen')

    room_dimension = [room_l, room_w, room_h];
    array_center = [array_x, array_y, array_z]; 
    source_position = [source_position_x, source_position_y, source_position_z];
    
    % eigenMike coordinates. The left-down corner is original point.
    MicCoor = MicCoorRel + array_center;       

    S = 2 * (room_h*room_w+room_w*room_l+room_l*room_h);
    V = room_h * room_w * room_l;
    alpha = 0.161 * V / (S * T60);                      
    beta = sqrt(1-alpha);

    % represent array center, room, sound source in samples
    array_center_in_sample = array_center /c * fs;
    source_position_in_sample = source_position / c * fs;
    room_dimension_in_sample = room_dimension / c * fs;   
                
    ht = sroom_hrtf(array_center_in_sample, ...
                    source_position_in_sample, ...
                    room_dimension_in_sample, ... 
                    0, ...     % roomId, used param for now
                    beta, ... % reflection coefficient, unused for now
                    tf_len, ...
                    fs, ...
                    c);
%     ht = ht';
%     ht = data_preprocess(ht);
    ht = sortrows(ht, 1);  % sort
    if order ~= -1
        ht = ht(find(ht(:, 5) <= 3), :);
    end
    % generate TF
    TF_Matrix = CalculateTF(fs, tf_len, c, MicCoor, ...  
                            array_center, TF_Lib, ...
                            dist_in_lib, angle_in_lib, ...
                            beta, ht);
                    
    
    parsave1([rir_path,'azi_', num2str(row+1), '_ele_', num2str(col+1), '_', num2str(pos_base+1), '.mat'], struct('TF', TF_Matrix, 'room', room_dimension, 'array_center', array_center, ... 
    'source_center_dist', source_center_dist, 'azi', azi / pi * 180, 'ele', ele / pi * 180, 'T60', T60, 'beta', beta));
    
end

disp('TF generated finished...')

if(PARA)
    delete(gcp('nocreate'))
end


function num = draw(range)
    if range(1) > range(2)
        error('error\n');
    end
    num = (range(2)-range(1))*rand+range(1);
end

function cellArr = assignCell(cellArr, arr)
    if length(cellArr) ~= length(arr)
        error('error\n');
    end
    for idx = 1:length(cellArr)
        cellArr{idx} = arr(idx);
    end
end
