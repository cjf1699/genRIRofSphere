function TF_single = sphere_hrtf(Fs,TF_len,theta,distance)
% This function will calcuate the ANECHOIC impulse response between a far-field sound source 
% located at the given azimuth 'theta' and distance 'distance'.
% Fs: the sample rate at which the impulse response would be calculated.
% TF_len: the length of the calculated impulse response
% Note: In order to support the RIR calculation which considers multiple image sources, the distance attenuation is set to that of 1m 
freq_resolution = Fs/TF_len;  % 
freq_axis = 0:freq_resolution:Fs-freq_resolution; % 
% 
sph_radius = 0.042;	%sph_radius
sound_speed = 344;  %sound_speed
threshold = 0.0000001; %threshold

mid_pos = floor(TF_len/2)+1;
azimuth = theta;
azi_num = length(azimuth);  % 1
azi_num_half = floor(azi_num/2); % 0
airDensity = 16.367;
% 
TF_single = zeros(TF_len,azi_num);  % just calculate TF of one angle
for azi_i=1:azi_num_half+1          % loop only once
    tf_tem = zeros(mid_pos,1);
    for freq_i = 2:mid_pos
       h_ff = airDensity/(4*pi*1)*exp(1i*2*pi*freq_axis(freq_i)*1/sound_speed);
       tf_tem(freq_i) = conj(sphere(sph_radius, distance, azimuth, freq_axis(freq_i), sound_speed,threshold)*h_ff);%;
    end
    tf_tem(1) = abs(tf_tem(2));% 
    if mod(TF_len,2) == 0
        tf_tem(end) = real(tf_tem(end));
        TF_single(:,azi_i) = ifft([tf_tem; conj(flipud(tf_tem(2:end-1)))]);
    else
        TF_single(:,azi_i) = ifft([tf_tem; conj(flipud(tf_tem(2:end)))]);
    end
end
% TF_single(:,azi_num_half+2:end) = fliplr(TF_single(:,2:azi_num_half));

end
