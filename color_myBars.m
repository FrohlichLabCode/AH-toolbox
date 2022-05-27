function color_myBars(CT, numDiv)
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects
for i=1:numDiv
    if i == 1 % Set the color of the first box to gray
        box = a(3*numDiv+1-i); color = [0,0,0]; set(box,'Color',color);
    else
        box = a(3*numDiv+1-i);   
        color = CT(i, :);
        set(box, 'Color', color);
    end
end
end
