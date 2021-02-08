function [ TF_Lib ] = CalculateTF_Lib(fs, tf_len, theta, distance)
%   This function calculates a series of impulse responses(IRs) when sound source is located at the given 
%   theta and distance, forming an IR library.
%   Note: The path delay and path attenuation of TFs in this library are  ALL for
%   distance==1m, while how the TFs look like depends on their distances
%   respectively.
    TF_Lib = cell(length(distance),length(theta));
    [m, n] = size(TF_Lib);
    
    fprintf('Start calculating the Lib!\n');
    
    for i = 1:m
        for j = 1:n
            TF_Lib{i,j}.distance = distance(i);
            TF_Lib{i,j}.theta = theta(j);
            fprintf('I am calculating TF for %fm, %f, %d TFs remaining...\n', TF_Lib{i,j}.distance, TF_Lib{i,j}.theta / pi * 180, m*n-((i-1)*n+j));
            TF_Lib{i,j}.data = sphere_hrtf(fs, tf_len, TF_Lib{i,j}.theta, TF_Lib{i,j}.distance);
        end
    end    
end

