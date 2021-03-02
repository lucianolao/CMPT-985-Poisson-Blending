function output = imblend(source, mask, target, transparent)
%Source, mask, and target are the same size (as long as you do not remove
%the call to fiximages.m). You may want to use a flag for whether or not to
%treat the source object as 'transparent' (e.g. taking the max gradient
%rather than the source gradient).

[imh, imw, nb] = size(source);

k = 0; % k = number of variables (pixels that will blend)
mapping = zeros(imh,imw);
for y = 1:imh
	for x = 1:imw
        if mask(y,x) ~= 0
            % mapping the different variables
            k = k+1;
            mapping(y,x) = k;
        end
	end
end

e = 1;
b = zeros(k, nb); % 3 columns for RGB pictures
maxNumberCoefficients = 5 * k;
global i j v position;
i = zeros(1,maxNumberCoefficients);
j = zeros(1,maxNumberCoefficients);
v = zeros(1,maxNumberCoefficients);
position = 0;

% Compute indices i, j, v for matrix A, as well as the vector b
for y = 1:imh
    for x = 1:imw
        if (mapping(y,x) ~= 0)
            for ch = 1:nb
                boundary_t = 0;
                if (x ~= 1) && (x ~= imw) && (y ~= 1) && (y ~= imh)
                    % derivative on both axis
                    boundary_t = assignIJV(e, mapping(y,x), 4, target(y,x,ch), boundary_t, ch);
                    boundary_t = assignIJV(e, mapping(y,x+1), -1, target(y,x+1,ch), boundary_t, ch);
                    boundary_t = assignIJV(e, mapping(y,x-1), -1, target(y,x-1,ch), boundary_t, ch);
                    boundary_t = assignIJV(e, mapping(y+1,x), -1, target(y+1,x,ch), boundary_t, ch);
                    boundary_t = assignIJV(e, mapping(y-1,x), -1, target(y-1,x,ch), boundary_t, ch);
                    if transparent
                        if abs(source(y,x,ch)-source(y,x+1,ch)) > abs(target(y,x,ch)-target(y,x+1,ch))
                            b(e,ch) = b(e,ch) + source(y,x,ch) - source(y,x+1,ch);
                        else
                            b(e,ch) = b(e,ch) + target(y,x,ch) - target(y,x+1,ch);
                        end
                        if abs(source(y,x,ch)-source(y,x-1,ch)) > abs(target(y,x,ch)-target(y,x-1,ch))
                            b(e,ch) = b(e,ch) + source(y,x,ch) - source(y,x-1,ch);
                        else
                            b(e,ch) = b(e,ch) + target(y,x,ch) - target(y,x-1,ch);
                        end
                        if abs(source(y,x,ch)-source(y+1,x,ch)) > abs(target(y,x,ch)-target(y+1,x,ch))
                            b(e,ch) = b(e,ch) + source(y,x,ch) - source(y+1,x,ch);
                        else
                            b(e,ch) = b(e,ch) + target(y,x,ch) - target(y+1,x,ch);
                        end
                        if abs(source(y,x,ch)-source(y-1,x,ch)) > abs(target(y,x,ch)-target(y-1,x,ch))
                            b(e,ch) = b(e,ch) + source(y,x,ch) - source(y-1,x,ch);
                        else
                            b(e,ch) = b(e,ch) + target(y,x,ch) - target(y-1,x,ch);
                        end
                    else
                        b(e,ch) = 4*source(y,x,ch) - source(y,x+1,ch) - source(y,x-1,ch) - source(y+1,x,ch) - source(y-1,x,ch);
                    end
                    b(e,ch) = b(e,ch) + boundary_t;
                elseif (x ~= 1) && (x ~= imw)
                    % derivative on x
                    boundary_t = assignIJV(e, mapping(y,x), 2, target(y,x,ch), boundary_t, ch);
                    boundary_t = assignIJV(e, mapping(y,x+1), -1, target(y,x+1,ch), boundary_t, ch);
                    boundary_t = assignIJV(e, mapping(y,x-1), -1, target(y,x-1,ch), boundary_t, ch);
                    if transparent
                        if abs(source(y,x,ch)-source(y,x+1,ch)) > abs(target(y,x,ch)-target(y,x+1,ch))
                            b(e,ch) = b(e,ch) + source(y,x,ch) - source(y,x+1,ch);
                        else
                            b(e,ch) = b(e,ch) + target(y,x,ch) - target(y,x+1,ch);
                        end
                        if abs(source(y,x,ch)-source(y,x-1,ch)) > abs(target(y,x,ch)-target(y,x-1,ch))
                            b(e,ch) = b(e,ch) + source(y,x,ch) - source(y,x-1,ch);
                        else
                            b(e,ch) = b(e,ch) + target(y,x,ch) - target(y,x-1,ch);
                        end
                    else
                        b(e,ch) = 2*source(y,x,ch) - source(y,x+1,ch) - source(y,x-1,ch);
                    end
                    b(e,ch) = b(e,ch) + boundary_t;
                elseif (y ~= 1) && (y ~= imh)
                    % derivative on y
                    boundary_t = assignIJV(e, mapping(y,x), 2, target(y,x,ch), boundary_t, ch);
                    boundary_t = assignIJV(e, mapping(y+1,x), -1, target(y+1,x,ch), boundary_t, ch);
                    boundary_t = assignIJV(e, mapping(y-1,x), -1, target(y-1,x,ch), boundary_t, ch);
                    if transparent
                        if abs(source(y,x,ch)-source(y+1,x,ch)) > abs(target(y,x,ch)-target(y+1,x,ch))
                            b(e,ch) = b(e,ch) + source(y,x,ch) - source(y+1,x,ch);
                        else
                            b(e,ch) = b(e,ch) + target(y,x,ch) - target(y+1,x,ch);
                        end
                        if abs(source(y,x,ch)-source(y-1,x,ch)) > abs(target(y,x,ch)-target(y-1,x,ch))
                            b(e,ch) = b(e,ch) + source(y,x,ch) - source(y-1,x,ch);
                        else
                            b(e,ch) = b(e,ch) + target(y,x,ch) - target(y-1,x,ch);
                        end
                    else
                        b(e,ch) = 2*source(y,x,ch) - source(y+1,x,ch) - source(y-1,x,ch);
                    end
                    b(e,ch) = b(e,ch) + boundary_t;
                else
                    % corner
                    % adding 1 equation so matrix A still size K,K
                    boundary_t = assignIJV(e, mapping(y,x), 1, target(y,x,ch), boundary_t, ch);
                    b(e,ch) = source(y,x,ch) + boundary_t;
                end
            end % end loop for channels
            e = e+1;
        end % end if of mask ~= 0
    end
end


% Create matrix A from the indices
i = i(1:position);
j = j(1:position);
v = v(1:position);
A = sparse(i,j,v);

% Solve A\b, and copy new pixels to their respective coordinates
output = zeros(imh, imw, nb);
for ch = 1:nb
    solution = A\b(:,ch);
    % error = sum(abs(A*solution-b(:,ch)));
    % disp(error);
    for y = 1:imh
        for x = 1:imw
            if (mapping(y,x) ~= 0)
                output(y,x,ch) = solution(mapping(y,x));
                % output(y,x,ch) = solution(ch,mapping(y,x));
            else
                % mask == 0, just copy pixels from target image
                output(y,x,ch) = target(y,x,ch);
            end
        end
    end
end    

%output = source .* mask + target .* ~mask;

end % end of function


% Function to keep track of the indices for matrix A
function [boundary_t] = assignIJV(row, mapping, value, target_pixel, b, channel)
    global i j v position;
    if (mapping ~= 0)
        if (channel == 1) % otherwise it doesnt need to create matrix A again
            position = position + 1;
            i(position) = row;
            j(position) = mapping;
            v(position) = value;
        end
        boundary_t = b;
    else
        % mask == 0, get values from the target image
        boundary_t = b + target_pixel;
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As explained on the web page, we solve for output by setting up a large
% system of equations, in matrix form, which specifies the desired value or
% gradient or Laplacian (e.g.
% http://en.wikipedia.org/wiki/Discrete_Laplace_operator)

% The comments here will walk you through a conceptually simple way to set
% up the image blending, although it is not necessarily the most efficient
% formulation. 

% We will set up a system of equations A * x = b, where A has as many rows
% and columns as there are pixels in our images. Thus, a 300x200 image will
% lead to A being 60000 x 60000. 'x' is our output image (a single color
% channel of it) stretched out as a vector. 'b' contains two types of known 
% values:
%  (1) For rows of A which correspond to pixels that are not under the
%      mask, b will simply contain the already known value from 'target' 
%      and the row of A will be a row of an identity matrix. Basically, 
%      this is our system of equations saying "do nothing for the pixels we 
%      already know".
%  (2) For rows of A which correspond to pixels under the mask, we will
%      specify that the gradient (actually the discrete Laplacian) in the
%      output should equal the gradient in 'source', according to the final
%      equation in the webpage:
%         4*x(i,j) - x(i-1, j) - x(i+1, j) - x(i, j-1) - x(i, j+1) = 
%         4*s(i,j) - s(i-1, j) - s(i+1, j) - s(i, j-1) - s(i, j+1)
%      The right hand side are measurements from the source image. The left
%      hand side relates different (mostly) unknown pixels in the output
%      image. At a high level, for these rows in our system of equations we
%      are saying "For this pixel, I don't know its value, but I know that
%      its value relative to its neighbors should be the same as it was in
%      the source image".

% commands you may find useful: 
%   speye - With the simplest formulation, most rows of 'A' will be the
%      same as an identity matrix. So one strategy is to start with a
%      sparse identity matrix from speye and then add the necessary
%      values. This will be somewhat slow.
%   sparse - if you want your code to run quickly, compute the values and
%      indices for the non-zero entries in A and then construct 'A' with a
%      single call to 'sparse'.
%      Matlab documentation on what's going on under the hood with a sparse
%      matrix: www.mathworks.com/help/pdf_doc/otherdocs/simax.pdf
%   reshape - convert x back to an image with a single call.
%   sub2ind and ind2sub - how to find correspondence between rows of A and
%      pixels in the image. It's faster if you simply do the conversion
%      yourself, though.
%   see also find, sort, diff, cat, and spy


