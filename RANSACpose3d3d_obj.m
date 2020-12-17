function [T,inlierPts,inlierObjs] = RANSACpose3d3d_obj(objectPts1,objectPts2,obj_stat)
T= eye(4);
iterations = 0;
N = length(objectPts1);
objKinds = max(obj_stat(:,:));


foundT = 0;
inlierObjs = zeros(objKinds,1);
inlierPts = zeros(N,1);
objThreshold = 0.5;
inlierThreshold = 0.5;
errorThreshold = .3;
testPtsN = 4;

while (iterations < 200000 && foundT==0 )
    inlierObjs = zeros(objKinds,1);
    inlierPts = zeros(N,1);
%     disp(iterations)
    iterations = iterations+1;    
    testIndex = randperm(N,testPtsN);
        
    T_test = (objectPts1(testIndex,:)\objectPts2(testIndex,:))';
    
    totalKind = 0;
    majorKind = 0;
    for i = 1:objKinds
        objIdex = find(obj_stat(:)==i);
        objectPts_1 = objectPts1(objIdex,1:4)';
        objectPts_2 = objectPts2(objIdex,1:4)';
        if (~isempty(objectPts1))
            totalKind = totalKind+1;
            error_Pts = vecnorm(T_test*objectPts_1-objectPts_2);
            [GoodptsErroridx, GoodptsError] = find(error_Pts<errorThreshold);
            
%             error_object = mean(GoodptsError);
            error_object = mean(error_Pts);
            if (error_object<errorThreshold || length(GoodptsError)/length(error_Pts) > inlierThreshold)
%             if (error_object<errorThreshold)
                inlierPts(objIdex(GoodptsErroridx==1)) = 1;
                majorKind = majorKind+0.5;
                inlierObjs(i) = 1;
                if (i==1) %% background weight is higher
                    majorKind = majorKind+0.5;
                end
%                 disp(T_test)
%                 disp(i)
            end
        end
    end
    
%     error_ttPts = vecnorm(T_test*objectPts1'-objectPts2');
%     [inlierPts, inlierPtsError]= find(error_ttPts<errorThreshold);
    inlierNum = length(find(inlierPts(:)>0));
    if (majorKind/totalKind>objThreshold && inlierNum>20 )
        foundT = 1;
        T = T_test;
        disp("iterations "+iterations)
        disp("inlier " + inlierNum)
        disp("inlier ratio " + inlierNum/N )
    end
    
%     disp(T_test)
%     disp(inlierObjs')
%     disp(inlierNum)

end

if (foundT == 0)
    disp("failed")
%     disp(T_test)
%     disp(inlierObjs')
%     disp(inlierNum)
end


end
