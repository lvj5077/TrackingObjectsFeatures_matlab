for i = 1:7
        objIdex = find(obj_stat(:)==i);
        objectPts_1 = vpts1(objIdex,:);
        objectPts_2 = vpts2(objIdex,:);
        
        if (~isempty(objectPts_1))
            k =1;

            error_Pts = sampsonErrf(F_test, objectPts_1, objectPts_2);
            [GoodptsErroridx, GoodptsError] = find(error_Pts<k*SamponDistThreshold);
            rr = length(GoodptsError)/length(error_Pts)
            if (rr > feature_inlierRatio)
     
%                 inlierPts(objIdex(GoodptsErroridx==1)) = 1;
%                 if (i==1) %% background weight is higher
%                     static_weight = static_weight+0.5*rr;
%                 end
%                 static_weight = static_weight+0.5*rr;
%                 inlierObjs(i) = 1;

            end
        end
        
       
end
    