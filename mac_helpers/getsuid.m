function suid = getsuid(varargin)

nostop = 0;
getArgs(varargin,{'nostop=0'});

suid = mglGetParam('sunetID');

if isempty(suid)
    warning('No suid was found in mglGetParam');
    if nostop
        suid = 'temp';
    else
        in = input('You can set a suid now or [enter] to use a temporary id: ','s');
        if ~isempty(in)
            mglSetParam(in);
            suid = in;
        else
            suid = 'temp';
        end
    end
end
