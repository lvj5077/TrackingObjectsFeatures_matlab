function [F,inlierPts,inlierObjs] = RANSACpose2d2d_obj(vpts1,vpts2,obj_stat)

% vpts1 = normalize(vpts1);
% vpts2 = normalize(vpts2);
F= zeros(3);
iterations = 0;
N = length(vpts1);
objKinds = max(obj_stat(:,:));

feature_inlierRatio=0.4;
hpNum=8;
confidence=0.99;
iterationMax=round(log(1-confidence)/log(1-(1-feature_inlierRatio)^hpNum));

iterationMax = max(iterationMax,5000);
foundT = 0;
inlierObjs = zeros(objKinds,1);
inlierPts = zeros(N,1);
% inlierPts = ones(N,1);
staticLevelThreshold = 0.8;
SamponDistThreshold = 0.3;

static_weight_localBest = 0;

while (iterations < iterationMax && foundT==0 )
    inlierObjs = zeros(objKinds,1);
    inlierPts = zeros(N,1);
    iterations = iterations+1;  
    
    testIndex = randperm(N,hpNum);
        
%     T_test = fundamental_matrix(x1, x2);
    [F_test,~,~] = estimateFundamentalMatrix(vpts1(testIndex,:),vpts2(testIndex,:),'Method','Norm8Point');
    
    all_weight = 0;
    static_weight = 0;

    for i = 1:objKinds
        objIdex = find(obj_stat(:)==i);
        objectPts_1 = vpts1(objIdex,:);
        objectPts_2 = vpts2(objIdex,:);
        
        if (~isempty(objectPts_1))
            k =1;
            all_weight = all_weight+0.5;
            if (i ==1)
                k=2;
                all_weight = all_weight+0.5;
            end
            
            error_Pts = sampsonErrf(F_test, objectPts_1, objectPts_2);
            [GoodptsErroridx, GoodptsError] = find(error_Pts<k*SamponDistThreshold);
            rr = length(GoodptsError)/length(error_Pts);
            if (rr > feature_inlierRatio)
                inlierPts(objIdex(GoodptsErroridx==1)) = 1;
                if (i==1) %% background weight is higher
                    static_weight = static_weight+0.5;
                end
                static_weight = static_weight+0.5;
                inlierObjs(i) = 1;

            end
        end
        
       
    end
    
    if (static_weight/all_weight > static_weight_localBest)
        static_weight_localBest = static_weight/all_weight;
        inlierPts_localBest = inlierPts;
        inlierObjs_localBest = inlierObjs;
        F_test_localBest = F_test;
    end

    if (static_weight/all_weight >staticLevelThreshold)
        foundT = 1;
        F = F_test;
        all_weight
        disp("iterations "+iterations)
%         disp("inlier " + inlierNum)
%         disp("inlier ratio " + inlierNum/N )
    end
    
%     disp(T_test)
%     disp(inlierObjs')
%     disp(inlierNum)

end

if (foundT == 0)
    all_weight
    disp("failed")
    disp(iterations)
    static_weight_localBest
    F = F_test_localBest;
    inlierPts = inlierPts_localBest;
    inlierObjs = inlierObjs_localBest;
%     disp(T_test)
%     disp(inlierObjs')
%     disp(inlierNum)
end


end
