function [ theta ] = CalculateTheta( centre, source, mic )
%   ���㾵��Դ���������ģ���˷�֮��ļн�
%   a: source, centre
%   b: mic, centre
%   c: source, mic
a = sqrt(sum((source - centre).^2, 2));  % a ����������Դ�������
b = sqrt(sum((mic - centre).^2, 2));
c = sqrt(sum((source - mic).^2, 2));
cos_theta = (a.^2+b^2-c.^2)./(2*a*b);
theta = acos(cos_theta);
end

