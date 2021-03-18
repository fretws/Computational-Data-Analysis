function trackPoints = trackByColor(frames,start,duration,colorMask)
%trackByColor Tracks the median location of the colors specified by the
%given colorMask function.
%   Returns a matrix of x,y coordinates of the median location of the
%   colors flagged by colorMask for every frame within frames from start to
%   start+duration-1, inclusive.
    trackPoints = zeros(duration, 2, 'uint16');
    for k=start:start+duration-1
        framek = frames(:,:,:,k);
        [y,x] = find(colorMask(framek));
        trackPoints(k-start+1, :) = [median(x), median(y)];
    end
end

