function TF_Matrix = CalculateTF(fs, tf_Len, sound_speed, MicCoor, ...
                    centre, TF_Lib, distances, angles, currBeta, ht)
 
% fs: sample rate
% tf_Len: the total points of impulse response
% sound speed
% MicCoor : absolute cartesian coordinates of Mics
% centre: absolute cartisian coordinate of array centre
% TF_Lib: pre-calculated TF library, to reduce computation
% distances: values of distances considered in the TF_Lib
% angles: values of angles considered in the TF_Lib
% currBeta: reflection coefficient
% ht: information of multiple image sources
    
    nMic = size(MicCoor, 1);
    
    numOfReflects = ht(:, 5);

    row = size(ht, 1);

    dist = sqrt(sum(ht(:,2:4).^2, 2));

    TF_perSourceLocResult = zeros(tf_Len, nMic); % given a source location, for each mic, summation over all  sources

    for mic_index = 1:nMic
%         fprintf('calculating mic %d...\n', mic_index);
%             CalculateTheta: check pass
        inner_angle = CalculateTheta(centre, ht(:,2:4)+centre, MicCoor(mic_index,:)); % 1 centre, many sources, 1 mic
        for source_index = 1:row                
%                 fprintf('HHHHH....calculating  source %d...\n', source_index);
            [disValue, disLocation] = min(abs(distances - dist(source_index)));
            [angValue, angLocation] = min(abs(angles - inner_angle(source_index)));
            
            TF_perSourceLoc = TF_Lib{disLocation, angLocation}.data * (currBeta.^numOfReflects(source_index)) / dist(source_index);
            
            % note: path delay  temp = [zeros(delay, 1)is needed to be compensated.
            if dist(source_index) < 1
                delay = round((1 - dist(source_index)) / sound_speed * fs);
            else
                delay = round(dist(source_index) / sound_speed * fs);
            end
            if delay > tf_Len
                continue
            end
            temp = TF_perSourceLoc;
            
            if dist(source_index) < 1
                temp = [temp(delay+1:end); zeros(delay, 1)];
            else
                temp = [zeros(delay, 1); temp(1:tf_Len - delay)];
            end
            TF_perSourceLoc = temp;
            
            TF_perSourceLocResult(:, mic_index)  = TF_perSourceLocResult(:, mic_index) + TF_perSourceLoc;
        end 
         
    end
    TF_Matrix = TF_perSourceLocResult;
    TF_Matrix = TF_Matrix / max(max(TF_Matrix));  % normalize
end
