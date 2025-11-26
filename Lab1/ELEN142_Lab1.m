x1 = [1,2,3,4];
x2 = [3,4,5,6] %no semicolon

x3 = [x1 x2];
x4 = [x1; x2];

length(x3)
size(x4)

x5 = 5*x4

%%
clear;
clc;
close all;

f = 5;
t = 1:0.01:2;
x = 3*cos(2*pi*f*t);

fy = 4;
ty = 1:0.01:2;
y = 2*sin(2*pi*fy*ty);

plot(t, x, t, y);
xlabel('Time');
ylabel('Amplitude');

%%
clear; clc; close all;
x = -10:0.01:10;

subplot(2,2,1)
y1 = sinc(x);
plot(x,y1)
title('Subplot 1: sinc(x)')

subplot(2,2,2)
y2 = sinc(x-2);
plot(x,y2)
title('Subplot 2: sinc(x-2)')

subplot(2,2,3)
sum = y1 + y2;
plot(x,sum)
title('Subplot 3: y1 + y2')

subplot(2,2,4)
difference = y1 - y2;
plot(x,difference)
title('Subplot 4: y1 - y2')

%%
clear; clc; close all;
sum = 0;
for index = 1:100
   sum = sum + index;
end
disp(sum);

%%
clear; clc; close all;
for c = 1:100
    if mod(c,5) == 0 && mod(c,3) == 0
        disp("FizzBuzz");
    elseif mod(c,5) == 0
        disp("Buzz");
    elseif mod(c,3) == 0
        disp("Fizz");
    else
        disp(c);
    end
end

%%
clear; clc; close all;

x = 1:1:100;
[evens,odds] = lab1_even_odd(x)


%%
clear; clc; close all;

Fs = 1000;
ts = 1/Fs;
t = 0:ts:10;
x1 = cos(2*pi*100*t);
x2 = cos(2*pi*200*t);
x = x1+x2;
X = fft(x);
shift = fftshift(X);
freqaxis = Fs*(linspace(-0.5,0.5,length(x)));
plot(freqaxis, abs(shift));
xlabel('Frequency');
ylabel('Amplitude');

%%
clear; clc; close all;

[y, Fs] = audioread('goat.wav');
sound(y,Fs); % listen to the audio

%%
clear; clc; close all;

[y, Fs] = audioread('goat.wav');
sound(y,Fs); % listen to the audio
Y = fft(y);
shift = fftshift(Y);
freqaxis = Fs*(linspace(-0.5,0.5,length(y)));
plot(freqaxis, abs(shift));
xlabel('Frequency');
ylabel('Amplitude');

%%
function [evens, odds] = lab1_even_odd(x)
    evens = x(mod(x,2)==0);
    odds = x(mod(x,2)==1);
end