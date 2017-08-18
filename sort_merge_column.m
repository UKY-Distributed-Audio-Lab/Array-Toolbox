function sorted = sort_merge_column( data, column )
%  Sorts an array according to the values in COLUMN, from least to greatest

lenD = size(data,1);
if(lenD<2)
  sorted = data;
else
  middle = cast(floor(lenD/2),'uint16');
  L = data(1:middle,:);
  R = data(middle+1:end,:);
  L = sort_merge_column(L,column);
  R = sort_merge_column(R,column);
  sorted = merge(L, R,column);
end

function dataOut = merge(L, R, column)
lenL = size(L,1);
lenR = size(R,1);
i = 0;
j = 0;
merged = [];
while(i<lenL||j<lenR)
  if (i<lenL && j<lenR)
    if(L(i+1,column)<=R(j+1,column))
      merged(i+j+1,:) = L(i+1,:);
      i = i+1;
    else
      merged(i+j+1,:) = R(j+1,:);
      j = j+1;
    end
  elseif(i<lenL)
    merged(i+j+1,:) = L(i+1,:);
    i = i+1;
  elseif(j<lenR)
    merged(i+j+1,:) = R(j+1,:);
    j = j+1;
  end
end
dataOut = merged;