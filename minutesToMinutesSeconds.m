function [mins,secs] = minutesToMinutesSeconds(minutes)
%MINUTESTOMINUTESSECONDS converts a non integer number of minutes to
%   minutes and seconds
%   Returns number of minutes and number of seconds
    mins = floor(minutes);
    secs = (minutes-mins)*60;
end

