%lowpass,weightmap�������²�������
function I=extDownsample(lowpass,filter)
I=lowpass;
for l = 1:2
    I = downsample(I,filter);
end

end