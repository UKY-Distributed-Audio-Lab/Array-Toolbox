function tdoaplot(source, tstmp, mpos)

% plots the detected source in a circular cage
% source is given as xy coords. timestamp TSTMP in sec. and mpos in xyz
% coords. w/ each column corresponding to a different column.

draw_circle(0,0, .2096, 'black');
scatter(0,0, 30, 'magenta');
scatter( mpos(1,:), mpos(2,:), 60, 'blue');
scatter( source(1,1),source(2,1),60, 'red');
