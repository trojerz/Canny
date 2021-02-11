function Canny_detector(record)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this implementation of Canny edge detector takes a record name as an
    % argument: Canny_detector '0001.png' 

    % set the parameter improve to 1 if you want improved version, otherwise 0
    improve = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load the image  and transform in the correct format
    CT_image = imread (record);    
    
    if improve
        % figure out that lower constant is better, but not too much
        % the best if somewhere from 0 to 0.2
        % const = rand(1) / 5
        % after some testing, our constant should be around 0.15
        %const = 0.15;
        
        const = 0.15;
        
        % for brains -> around 0.1
        image_mean_value = mean2(CT_image);
        image_sd = std2(CT_image);
        lower_bound = (image_mean_value - image_sd);
        upper_bound = (image_mean_value + image_sd);
        threshold_l = max(1,const * lower_bound) / 255;
        threshold_h = min(254,const * upper_bound) / 255;
    else
        % this threshold was determined empirically 
        threshold_l = 0.060;
        threshold_h = 0.160;
    end
    
    CT_image = im2double(CT_image);
    matlab_canny = edge(CT_image, 'Canny');
    %threshold_l = low;
    %threshold_h = high;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gaussian filter
    
    
    
    % like in the literature
    % 
    % const = 1/159;
    % B = const.*[2, 4, 5, 4, 2; 
    %             4, 9, 12, 9, 4;
    %             5, 12, 15, 12, 5;
    %             4, 9, 12, 9, 4;
    %             2, 4, 5, 4, 2 ];
    % 
    sigma = min(size(CT_image)) * 0.005;
    B = calcGauss(sigma);
    % need same size as CT_image - apply the filter
    CT_image_gauss = conv2(CT_image, B, 'same');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sobel kernel
    Bx = [-1,0,1;
          -2,0,2;
          -1,0,1];
    By = [1, 2, 1;
          0, 0, 0;
         -1,-2,-1];
    % applying sobel kernel
    CT_sobel_x = conv2(CT_image_gauss, Bx, 'same');
    CT_sobel_y = conv2(CT_image_gauss, By, 'same');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % edge graditent
    edge_gradient = sqrt((CT_sobel_x.^2) + (CT_sobel_y.^2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % angle
    angle = atan2(CT_sobel_y, CT_sobel_x );
    const_2 = 180 / pi;
    % angle in degrees
    angle_d = const_2 * angle;
    [size_x, size_y] = size(CT_image_gauss);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % all endge direction angle must be rounded to one of four angles
    % representing vertical, horizontal and two diagonals
    % 0, 45, 90, 135
    % negative angles -> positive angles with function make_positive
    for i = 1  : size_x
        for j = 1 : size_y 
            if ((make_positive(angle_d(i, j)) >= 0 ) && (make_positive(angle_d(i, j)) < 22.5) || (make_positive(angle_d(i, j)) >= 157.5) && (make_positive(angle_d(i, j)) < 202.5) || (make_positive(angle_d(i, j)) >= 337.5) && (make_positive(angle_d(i, j)) <= 360))
                angle_d(i, j) = 0;
            elseif ((make_positive(angle_d(i, j)) >= 22.5) && (make_positive(angle_d(i, j)) < 67.5) || (make_positive(angle_d(i, j)) >= 202.5) && (make_positive(angle_d(i, j)) < 247.5))
                angle_d(i, j) = 45;
            elseif ((make_positive(angle_d(i, j)) >= 67.5 && make_positive(angle_d(i, j)) < 112.5) || (make_positive(angle_d(i, j)) >= 247.5 && make_positive(angle_d(i, j)) < 292.5))
                angle_d(i, j) = 90;
            elseif ((make_positive(angle_d(i, j)) >= 112.5 && make_positive(angle_d(i, j)) < 157.5) || (make_positive(angle_d(i, j)) >= 292.5 && make_positive(angle_d(i, j)) < 337.5))
                angle_d(i, j) = 135;
            end;
        end;
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % non-maximum suppression
    % every pixel needs to be checked if it is a local maximum in the gradient
    % direction
    % check pixel by pixel for the gradient value, compared to neighbours in 
    % two directions. For example, if we are on the 'horizonal line', we check
    % neighbours North and South and compare the values. If current is greater,
    % we know that this is the edge!
    edges = zeros(size_x, size_y);
    % calculation for for each pixel 
    for i = 2 : size_x - 1 %careful for edges!
        for j = 2 : size_y - 1
            current = edge_gradient(i,j);
            if (angle_d(i,j) == 0) % horizontal line
                % check if current has greater gradient than one N / S pixel
                edges(i,j) = (current == max([current, edge_gradient(i,j+1), edge_gradient(i,j-1)])) * current;
            elseif (angle_d(i,j) == 45) % positive diagonal
                 % check if current has greater gradient than one NW / SE pixel
                edges(i,j) = (current == max([current, edge_gradient(i+1,j-1), edge_gradient(i-1,j+1)])) * current;
            elseif (angle_d(i,j) == 90) % vertical line
                 % check if current has greater gradient than one W / E pixel
                edges(i,j) = (current == max([current, edge_gradient(i+1,j), edge_gradient(i-1,j)])) * current;
            elseif (angle_d(i,j) == 135) % negative diagonal
                  % check if current has greater gradient than one SW / NE pixel
                edges(i,j) = (current == max([current, edge_gradient(i+1,j+1), edge_gradient(i-1,j-1)])) * current;
            end;
        end;
    end;
    % adjust the threshold to the image
    threshold_l = threshold_l * max(max(edges));
    threshold_h = threshold_h *  max(max(edges));

    final_CT = zeros(size_x, size_y);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Double thresholding & hysteresis -> connectiong edges
    for i = 1  : size_x
        for j = 1 : size_y        
            if (edges(i, j) > threshold_h)
                final_CT(i, j) = 255; %strong edge pixel
            elseif (edges(i, j) < threshold_l)
                final_CT(i, j) = 0; % suppressed
            else
                % conditions for weak edge, compare neighbour positions
                % 8 connected components. If at least one true, then weak edge is
                % strong edge
                east_pos = (edges(i+1,j) > threshold_h);
                west_pos = (edges(i-1,j) > threshold_h);
                north_pos = (edges(i,j+1) > threshold_h);
                south_pos = (edges(i,j-1) > threshold_h);
                sw_pos = (edges(i-1, j-1) > threshold_h);
                nw_pos = (edges(i-1, j+1) > threshold_h);
                ne_pos = (edges(i+1, j+1) > threshold_h);
                se_pos = (edges(i+1, j-1) > threshold_h);
                if (west_pos || east_pos || north_pos || south_pos || sw_pos || nw_pos || ne_pos || se_pos)
                    % weak edge becomes strong edge, if connected to strong
                    % edge
                    final_CT(i,j) = 255;
                end;
            end;
        end;
    end;
    % black - background
    % white - 'sceleton'
    final_CT = uint8(final_CT);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write the image
    imwrite(final_CT, strcat('./results/', record));
    
    %imwrite(matlab_canny, strcat('./results/', record));
    
    matlab_canny = matlab_canny .* 255;
    matlab_canny = uint8(matlab_canny);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % showing images
    % original
    figure;
    imshow(CT_image);
    colormap gray;
    % after detecting edges
    figure;
    imshow(edges);
    colormap gray;
    % after edge linking
    figure;
    imshow(final_CT);
    % matlab implementation
    %figure;
    %imshow(matlab_canny);