function y_obj = objectives(x_var)
global numVar AdjacentMatrix node

% copygene si : to store the k clusters and the elements
si=cell(0);
copygene=x_var;
for i=1:numVar
    %//store every component in each cluster
    if copygene(i) ~= -1
        s=[];
        s=[s i];%s.push_back(i);
        for j=i+1:numVar
            if copygene(i) == copygene(j)
                s=[s j];%s.push_back(j);
                copygene(j) = -1;
            end
        end % end j
        copygene(i) = -1;
        si=[si s];
    end %end if
end %end i
clusters = size(si,2);
% if SignedFlag == 0
    Temp_RA = 0.0; 
    Temp_RC = 0.0;
    for i =1:clusters
        vs_i = 0;
		ki_out = 0;
        siisize=size(si{i},2);
        for j = 1:siisize
            a = si{i}(j);
			kj_in = 0;
            for k=1:siisize
                b = si{i}(k);
				kj_in = kj_in + AdjacentMatrix(a,b);
            end %end k
            vs_i = vs_i + kj_in;
			ki_out = ki_out + (node(a).degree - kj_in);    
        end %end j   
        Temp_RA = Temp_RA + 1.0 * vs_i/siisize;
		Temp_RC = Temp_RC + 1.0 * ki_out/siisize;	
    end %end i
%     copygene=[];
%     si={};
    % /* for minimizing the objectives */
	y_obj(1) = 2 * (numVar - clusters) - Temp_RA;	% //KKM
	y_obj(2) = Temp_RC;	
% else if SignedFlag == 1 

end
                
                
                
                
                
                
                
                
                
                
                
                