function [V, S] = makeBox(bboxL, bboxW, scale)
    center = [bboxL/2, bboxW/2];
    dimension = [bboxL/scale, bboxW/scale];
    TR = center + dimension/2;
    BL = center - dimension/2;
    V = [BL;TR(1) BL(2);  TR; BL(1) TR(2);  ];
    S = [1 2; 2 3; 3 4; 4 1];
end

