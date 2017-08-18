function join_vars(varlabel, varstr, outfile)
%  This function will join the values of a same-named variable from as many
%  .mat files as the user selects. If user wants the variable to be sorted,
%  the user must specify a column to sort by. The result will be saved in
%  OUTFILE.

flag='YES';
varhold={};
%varname= genvarname(var);

while flag(1)=='Y'
   [fn, pn] = uigetfile('*.mat','Select .mat file.');
   matfile = [pn, fn];
   load(matfile);
   varhold= [varhold, varlabel];
   flag = questdlg('Would you like to select more files', 'Grpahics Flag', 'Yes' ,'No', 'No');
end

varFinal= [];
for i=1:size(varhold,2)
    varFinal= [varFinal; varhold{1,i}];
end

column_sort= 1;
varFinal= sort_merge_column(varFinal, column_sort);
save( outfile, 'varFinal' );
