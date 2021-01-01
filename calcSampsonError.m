function d = calcSampsonError(F, x1, x2)
% x1 = normalise(x1);
% x2 = normalise(x2);  
% normalise epipolar lines (unit normal)

l1 = (F' * x2);
l2 = (F * x1);
l1 = l1 ./ sqrt(repmat(l1(1,:).^2 + l1(2,:).^2 + l2(1,:).^2 + l2(2,:).^2, size(l1,1), 1));
l2 = l2 ./ sqrt(repmat(l1(1,:).^2 + l1(2,:).^2 + l2(1,:).^2 + l2(2,:).^2, size(l1,1), 1));

% distance from points x1 to epipolar line of corresponding x2 in I1
d1 = abs(sum(l1.*x1));
% distance from points x2 to epipolar line of corresponding x1 in I2
d2 = abs(sum(l2.*x2));
% sum of distances
d = d1 + d2; 
end