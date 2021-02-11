function angle_positive = make_positive(angle)
    if angle < 0
        angle_positive = angle + 360;
    else
        angle_positive = angle;
    end
end