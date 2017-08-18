function srpvarsanalyze( output )
% This function takes in the output from srpvarstest.m and creates a
% structure of only 'true' detections. It then creates plots to show how
% different values of a variable correlate to detection rates. 

output_det=[];
titles= { 'fsb' 'nwlen' 'fc' 'fsresample' 'beta' };

for i=1:size(output,1)
    if output(i,7)==1
        output_det=[output_det; output(i,:)];
    end
end

for i=2:6
    %analyze for changes in window length nwlen
    dat= [unique(output_det(:,i)), zeros(size(unique(output_det(:,i)),1),1)];
    %iterate through output_det and increment dat(n,2) appropriately
    for j=1:size(output_det,1)
        n=find(output_det(j,i)==dat(:,1));
        dat(n,2)= dat(n,2)+1;
    end
    figure; 
    scatter(dat(:,1),dat(:,2));
    lims=ylim;
    ylim([0 lims(2)]);
    title( ['# of detections as ' titles{i-1} ' changes:'] );
end
end

