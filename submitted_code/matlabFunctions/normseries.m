function outputs = normseries(x)
 meanv = mean(x,2);
 tempx = x';
 varv = var(double(tempx));
 meanV = repmat(meanv,1,size(x,2));
 varV = repmat(sqrt(varv'), 1, size(x,2));
 
 meancorrected = double(x) - double(meanV);
 temp_outputs = meancorrected./varV;
 temp_outputs(isnan(temp_outputs)) = 0;
 outputs = temp_outputs;
 end
 