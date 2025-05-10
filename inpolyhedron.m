function in = inpolyhedron(varargin)

    narginchk(2,3);

    if nargin == 2
        if ~isstruct(varargin{1}) || ~all(isfield(varargin{1},{'vertices','faces'}))
            error('First argument must be a structure with "vertices" and "faces" fields, or call with F, V, P.');
        end
        fv = varargin{1};
        node = fv.vertices;
        face = fv.faces;
        testp = varargin{2};
    else
        face = varargin{1};
        node = varargin{2};
        testp = varargin{3};
    end

    if size(node,2) ~= 3
        error('Vertices must be an Nv-by-3 matrix.');
    end
    if size(face,2) ~= 3
        error('Faces must be an Nf-by-3 matrix (triangles).');
    end
    if size(testp,2) ~= 3
        error('Test points must be an Np-by-3 matrix.');
    end

    num_points = size(testp, 1);
    in = false(num_points, 1);

    for p_idx = 1:num_points
        pt = testp(p_idx, :);
        intersections = 0;


        for f_idx = 1:size(face, 1)
            v1 = node(face(f_idx, 1), :);
            v2 = node(face(f_idx, 2), :);
            v3 = node(face(f_idx, 3), :);

           
            y1 = v1(2) - pt(2); z1 = v1(3) - pt(3);
            y2 = v2(2) - pt(2); z2 = v2(3) - pt(3);
            y3 = v3(2) - pt(2); z3 = v3(3) - pt(3);

        
            if ((y1 > 0 && y2 <= 0) || (y1 <= 0 && y2 > 0)) && ...
               (z1 + (y1 / (y1 - y2)) * (z2 - z1) > 0)
                intersections = intersections + 1;
            end
       
            if ((y2 > 0 && y3 <= 0) || (y2 <= 0 && y3 > 0)) && ...
               (z2 + (y2 / (y2 - y3)) * (z3 - z2) > 0)
                intersections = intersections + 1;
            end
           
            if ((y3 > 0 && y1 <= 0) || (y3 <= 0 && y1 > 0)) && ...
               (z3 + (y3 / (y3 - y1)) * (z1 - z3) > 0)
                intersections = intersections + 1;
            end
        end
        
       
        if mod(intersections, 2) == 1
            in(p_idx) = true;
        end
    end
    
 end % function inpolyhedron
