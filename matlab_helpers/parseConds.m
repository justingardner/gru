
function pull = parseConds(pull)
% by Dan (2016)
%
%   A generic function for parsing the volume/condition labels that are
%   turned from getStimvol. Parses the "condition" cell to create an array
%   that includes the actual values (they have to be numeric).
%

pfields = fields(pull);

for pi = 1:length(pfields)
    dat = pull.(pfields{pi});
    
    val = [];
    for i = 1:length(dat.conds)
        name = dat.conds{i};
        val(end+1) = str2num(name(strfind(name,dat.str)+length(dat.str)+1:end));
    end
    dat.val = val;
    
    pull.(pfields{pi}) = dat;
end