%
temp = [-10:10:40]; % temperature values
humidity = [5:20:95]; % Humidity values
p = 38;  % Pressure values
f = 10e3;  % frequency in Hz
dist = 2;
for k=1:length(temp)
    for n=1:length(humidity)
        c1(k,n)= SpeedOfSound(temp(k),humidity(n),p);
        a1(k,n) = atmAtten(temp(k),p,humidity(n),dist,f);
    end
end