mglSetParam('fwAPIKey', 'cni.flywheel.io:2trpsDM6LhQq1Vv6cD');
mglSetParam('sunetID', 'jlg');
dbstop if error
%javaaddpath('/Users/justin/proj/gru/flywheel-sdk/api/rest-client.jar');
disp('hello world');
mrSetPref('maxBlockSize',4000000000);